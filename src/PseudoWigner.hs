{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts#-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE EmptyDataDecls #-}

module PseudoWigner(WindowFunc(..), makeWindow, pWignerVille) where

import Hilbert
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

data WindowFunc = Rect | Sin | Lanczos | Hanning | Hamming
  deriving (Read, Show)

makeWindow :: (A.RealFloat e, A.FromIntegral Int e) => WindowFunc -> A.Acc (Scalar Int) -> A.Acc (A.Array A.DIM1 e)
makeWindow func leng = 
  case func of 
    Rect -> A.fill (A.index1 $ A.the leng) 1.0
    Sin  -> A.generate (A.index1 $ A.the leng) (\sh -> let (A.Z A.:.x) = A.unlift sh in sin (pi*(A.fromIntegral x)/(A.fromIntegral $ A.the leng - 1)))
    Lanczos -> A.generate (A.index1 $ A.the leng) (\sh -> let (A.Z A.:.x) = A.unlift sh in sinc ((2*(A.fromIntegral x)/(A.fromIntegral $ A.the leng - 1)) - 1.0)) 
    Hanning -> A.generate (A.index1 $ A.the leng) (\sh -> let (A.Z A.:.x) = A.unlift sh in 0.5 - (0.5 * (cos (2*pi*(A.fromIntegral (x + 1))/(A.fromIntegral $ A.the leng + 1)))))
    Hamming -> A.generate (A.index1 $ A.the leng) (\sh -> let (A.Z A.:.x) = A.unlift sh in 0.54 - (0.46 * (cos (2*pi*(A.fromIntegral (x + 1))/(A.fromIntegral $ A.the leng + 1)))))

pWignerVille :: (A.Acc (A.Array A.DIM1 Float), A.Acc (A.Array A.DIM1 (ADC.Complex Float))) -> A.Acc (A.Array A.DIM2 Float)
pWignerVille (window, arr) = 
  let times = A.enumFromN (A.index1 leng) 0 :: A.Acc (Array DIM1 Int)
      leng = A.length arr
      taumx = taumaxs times window
      lims = limits taumx
  in A.map ADC.real $ A.transpose $ AMF.fft1D_m' AMF.Forward (A.Z A.:. 16 A.:. 16) $ createMatrix arr window taumx lims

taumax :: A.Exp Int -> A.Exp Int -> A.Exp Int -> A.Exp Int
taumax leng lh t = min (min (min t (leng - t - 1) ) (A.round (((A.fromIntegral leng :: A.Exp Double)/2.0) - 1))) lh

taumaxs :: A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Float) -> A.Acc (A.Array A.DIM1 Int)
taumaxs times window = 
  let leng = A.length times
      lh = (A.length window - 1) `div` 2
  in A.map (taumax leng lh) times                  

times :: Elt a => A.Acc (A.Array A.DIM1 a) -> A.Acc (A.Array A.DIM1 Int)
times arr = 
  let leng = A.length arr 
  in A.enumFromN (A.index1 leng) 0 :: A.Acc (Array DIM1 Int)

limits :: A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int)
limits taumaxs = 
  let funk = (\x -> 2*x + 1)
  in A.map funk taumaxs

moveUp ::  A.Acc (A.Array A.DIM1 Int) -> A.Exp Int -> A.Exp DIM2 -> A.Exp DIM2
moveUp taumaxs leng sh = 
  let taum t = taumaxs A.!! t 
  in (\(x,t) -> A.index2 ((x+(taum t)) `A.mod` leng) t) $ A.unlift $ A.unindex2 sh

generateValue :: A.Acc (A.Array A.DIM1 (ADC.Complex Float)) -> A.Exp Int -> A.Exp Int -> A.Exp Float -> A.Exp (ADC.Complex Float)
generateValue arr time tau h = (makeComplex h) * (arr A.!! (time + tau)) * (ADC.conjugate $ arr A.!! (time - tau))


createMatrix :: A.Acc (A.Array A.DIM1 (ADC.Complex Float)) -> A.Acc (A.Array A.DIM1 Float) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM2 (ADC.Complex Float)) 
createMatrix arr window taumaxs lims = A.transpose $ A.backpermute (A.index2 leng leng) (moveUp taumaxs leng) raw 
  where
    raw = A.generate (A.index2 leng leng) (\sh -> let (A.Z A.:.x A.:. t) = A.unlift sh
                                                      lim = lims A.!! t
                                                      taum = taumaxs A.!! t
                                                      h = window A.!! (lh + (x - taum))
                                                  in gen x t lim taum h)
    leng = A.length arr
    lh = (A.length window - 1) `div` 2 
    gen x t lim taum h = A.cond (x A.< lim) (generateValue arr t (x - taum) h) 0

epsilon :: RealFloat a => a
epsilon = encodeFloat 1 (fromIntegral $ 1-floatDigits epsilon)

sinc :: (A.RealFloat a) => A.Exp a -> A.Exp  a
sinc x =
   if abs x >= taylor_n_bound
     then sin x / x
     else 1 - x^2/6 + x^4/120
 where
  taylor_n_bound = sqrt $ sqrt epsilon