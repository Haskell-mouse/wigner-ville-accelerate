{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts#-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE EmptyDataDecls #-}

module Wigner where

--import RwData
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



wignerVille :: A.Acc (A.Array A.DIM1 (ADC.Complex Float)) -> A.Acc (A.Array A.DIM2 Float)
wignerVille arr = 
  let times = A.enumFromN (A.index1 leng) 0 :: A.Acc (Array DIM1 Int)
      leng = A.length arr 
      taumx = taumaxs times
      lims = limits taumx
  in A.map ADC.real $ A.transpose $ AMF.fft1D_m' AMF.Forward (A.Z A.:. 16 A.:. 16) $ createMatrix arr taumx lims 

taumax :: A.Exp Int -> A.Exp Int -> A.Exp Int
taumax leng t = min (min t (leng - t - 1) ) (A.round (((A.fromIntegral leng)/2.0) - 1.0 :: A.Exp Float))

taumaxs :: A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int)
taumaxs times = 
  let leng = A.length times
  in A.map (taumax leng) times                  

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

generateValue :: A.Acc (A.Array A.DIM1 (ADC.Complex Float)) -> A.Exp Int -> A.Exp Int -> A.Exp (ADC.Complex Float)
generateValue arr time tau = (arr A.!! (time + tau)) * (ADC.conjugate $ arr A.!! (time - tau))


createMatrix :: A.Acc (A.Array A.DIM1 (ADC.Complex Float)) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM2 (ADC.Complex Float)) 
createMatrix arr taumaxs lims = A.transpose $ A.backpermute (A.index2 leng leng) (moveUp taumaxs leng) raw 
  where
    raw = A.generate (A.index2 leng leng) (\sh -> let (A.Z A.:.x A.:. t) = A.unlift sh
                                                      lim = lims A.!! t
                                                      taum = taumaxs A.!! t
                                                  in gen x t lim taum)
    leng = A.length arr
    gen x t lim taum = A.cond (x A.< lim) (generateValue arr t (x - taum)) 0