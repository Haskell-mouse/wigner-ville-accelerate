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
--import qualified Data.Array.Accelerate.Interpreter as ALI
import qualified Data.Array.Accelerate.LLVM.Native as ALN
import qualified Data.Array.Accelerate.LLVM.PTX as ALP
import qualified Data.Array.Accelerate.Math.FFT as AMF
import qualified Data.Array.Accelerate.Data.Complex as ADC
import qualified Data.Vector.Storable as VS
import qualified Data.Array.Accelerate.IO as AI
import qualified Data.Vector as V 

{-
wignerVille :: A.Array A.DIM1 (ADC.Complex Double) -> V.Vector (A.Array A.DIM1 Double)
wignerVille arr = 
  let vectors = createRVectors arr
      dim1 = (A.Z A.:. (S.size $ shape arr))
  in V.map (ALN.run . (A.map ADC.real) . AMF.fft1D' AMF.Forward dim1) vectors -- `VST.using` (VST.parVector 1) -}

{-# INLINE taumax #-}
taumax :: Int -> Int -> Int
taumax leng t = min (min (t - 1) (leng - t) ) ((leng `mround` 2) -1)

{-# INLINE tau #-}
tau :: Int -> VS.Vector Int
tau taumax = let leng = (2*taumax + 1)
             in VS.force $ VS.generate leng (\x -> (-taumax + x))

{-# INLINE indices #-}
indices :: Int -> VS.Vector Int -> VS.Vector Int
indices leng tau = VS.force $ VS.map (((flip mod) leng) . ((+) leng)) tau


mround :: Int -> Int -> Int 
mround a !b = let half = b `div` 2
             in if (a `mod` b) >= half
                then ((a `div` b) + 1)
                else (a `div` b)

create :: A.Acc (A.Scalar Int) -> A.Acc (A.Array DIM2 Float)
create ss = 
  let sss = A.the ss
      sh = A.lift (A.Z A.:. sss A.:. sss)
  in A.generate sh (\ix -> let (A.Z A.:. (x :: A.Exp Int) A.:. y) = A.unlift ix in A.fromIntegral (y*10 + x))

makeFFT2 :: forall e. AMF.FFTElt e => A.Acc (A.Array A.DIM2 (ADC.Complex e)) -> A.Acc (A.Array A.DIM2 e)
makeFFT2 arr = A.map ADC.real $ AMF.fft1D_m' AMF.Forward (A.Z A.:. 16 A.:. 16) arr

type DIM2_INT = A.Z A.:. A.Exp Int A.:. A.Exp Int

makeGlobal ss = makeFFT2 $ makeComplex $ create ss