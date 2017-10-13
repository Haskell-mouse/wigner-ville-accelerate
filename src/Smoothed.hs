{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts#-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE EmptyDataDecls #-}

module Smoothed where --(WindowFunc(..), spWignerVille, matProduct, mkSmMatrix, taumaxs) where

import Hilbert
import qualified PseudoWigner as P
import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Array.Sugar as S 
import qualified Data.Array.Accelerate.Data.Fold as AF
import qualified Data.Array.Accelerate.Data.Monoid as AM
import qualified Data.Array.Accelerate.Interpreter as ALI
import qualified Data.Array.Accelerate.LLVM.Native as ALN
import qualified Data.Array.Accelerate.LLVM.PTX as ALP
import qualified Data.Array.Accelerate.Math.FFT as AMF
import qualified Data.Array.Accelerate.Data.Complex as ADC
import qualified Data.Vector.Storable as VS
import qualified Data.Array.Accelerate.IO as AI
import qualified Data.Vector as V 
import Debug.Trace

spWignerVille :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e) => 
  (A.Acc (A.Array A.DIM1 e), A.Acc (A.Array A.DIM1 (ADC.Complex e))) -> A.Acc (A.Array A.DIM2 e)
spWignerVille (window, arr) = 
  let times = A.enumFromN (A.index1 leng) 0 :: A.Acc (Array DIM1 Int)
      leng = A.length arr
      taumx = taumaxs times window
      lims = limits taumx
  in A.map ADC.real $ A.transpose $ AMF.fft1D_m' AMF.Forward (A.Z A.:. 16 A.:. 16) $ createFirstMatrix arr window taumx lims

taumax :: A.Exp Int -> A.Exp Int -> A.Exp Int -> A.Exp Int
taumax leng lh t = (min (min (min t (leng - t - 1) ) (A.round (((A.fromIntegral leng :: A.Exp Double)/2.0) - 1))) lh)

taumaxs :: (A.RealFloat e, Elt e) => 
  A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 e) -> A.Acc (A.Array A.DIM1 Int)
taumaxs times window = 
  let leng = A.length times
      lh = (A.length window - 1) `div` 2
  in A.map (taumax leng lh) times

taumaxs_ext :: A.Acc (A.Array A.DIM1 Int) -> A.Acc (Scalar Int) -> A.Acc (A.Array A.DIM1 Int)
taumaxs_ext taumaxs lG = A.map (\x -> x + A.the lG) taumaxs

times :: Elt a => A.Acc (A.Array A.DIM1 a) -> A.Acc (A.Array A.DIM1 Int)
times arr = 
  let leng = A.length arr 
  in A.enumFromN (A.index1 leng) 0 :: A.Acc (Array DIM1 Int)

limits :: A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int)
limits taumaxs = 
  let funk = (\taumax -> 2*taumax + 1)
  in A.map funk taumaxs

moveUp ::  A.Acc (A.Array A.DIM1 Int) -> A.Exp Int -> A.Exp DIM2 -> A.Exp DIM2
moveUp taumaxs leng sh = 
  let taum t = taumaxs A.!! t 
  in (\(x,t) -> A.index2 ((x+(taum t)) `A.mod` leng) t) $ A.unlift $ A.unindex2 sh

generateValue :: (A.RealFloat e, Elt e) => 
  A.Acc (A.Array A.DIM1 (ADC.Complex e)) -> A.Exp Int -> A.Exp Int -> A.Exp e -> A.Exp (ADC.Complex e)
generateValue arr time tau h = (arr A.!! (time + tau)) * (ADC.conjugate $ arr A.!! (time - tau))


createFirstMatrix :: (A.RealFloat e, Elt e) => 
  A.Acc (A.Array A.DIM1 (ADC.Complex e)) -> A.Acc (A.Array A.DIM1 e) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM2 (ADC.Complex e)) 
createFirstMatrix arr window taumaxs lims = A.backpermute (A.index2 leng leng) (moveUp taumaxs leng) raw 
  where
    raw = A.generate (A.index2 leng leng) (\sh -> let (A.Z A.:.x A.:. t) = A.unlift sh
                                                      lim = lims A.!! t
                                                      taum = taumaxs A.!! t
                                                      h = window A.!! (lh + (x - taum))
                                                  in gen x t lim taum h)
    leng = A.length arr
    lh = (A.length window - 1) `div` 2 
    gen x t lim taum h = A.cond (x A.< lim) (generateValue arr t (x - taum) h) 0

mkSmMatrix :: (A.RealFloat e, Elt e) => 
  A.Acc (A.Scalar Int) -> A.Acc (A.Array A.DIM1 e) -> A.Acc (A.Array A.DIM2 e)
mkSmMatrix arrLeng tWindow = 
  A.generate (A.index2 leng leng) (\sh -> let (A.Z A.:. x A.:. t) = A.unlift sh in gen x t)
  where 
    leng = A.the arrLeng
    wleng = A.length tWindow
    lG = (wleng - 1) `A.div` 2
    gen x t = A.cond  ((x A.> 2*t A.&& t A.<= lG) A.|| 
                       (x A.< t - (leng - t - 1) A.&& t A.>= leng - lG) A.|| 
                       (((x A.< t - lG) A.|| (x A.> t + lG)) A.&& (t A.> lG) A.&& (t A.< leng - lG))) 0 (tWindow A.!! (x + (wleng `A.div` 2) - t))

matProduct :: (A.RealFloat e, Elt e) => 
  A.Acc (A.Array A.DIM2 e) -> A.Acc (A.Array A.DIM2 e) -> A.Acc (A.Array A.DIM2 e)
matProduct wign smooth = 
  let A.Z A.:. wx A.:. _ = A.unlift $ A.shape wign :: A.Z A.:. A.Exp Int :. A.Exp Int 
      A.Z A.:. _ A.:. sy = A.unlift $ A.shape smooth :: A.Z A.:. A.Exp Int :. A.Exp Int
      aRep = A.replicate (A.lift $ (A.Z A.:. A.All A.:. sy A.:. A.All)) wign
      bRep = A.replicate (A.lift $ (A.Z A.:. wx A.:. A.All A.:. A.All)) (A.transpose smooth)
  in A.fold (+) 0 $ A.zipWith (*) aRep bRep 