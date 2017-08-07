{-# LANGUAGE BangPatterns #-}

module Wigner where

import RwData
import Hilbert
import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Array.Sugar as S 
import qualified Data.Array.Accelerate.Data.Fold as AF
import qualified Data.Array.Accelerate.Data.Monoid as AM
--import qualified Data.Array.Accelerate.Interpreter as ALI
import qualified Data.Array.Accelerate.LLVM.Native as ALN
--import qualified Data.Array.Accelerate.LLVM.PTX as ALP
import qualified Data.Array.Accelerate.Math.FFT as AMF
import qualified Data.Array.Accelerate.Data.Complex as ADC
import qualified Data.Vector.Storable as VS
import qualified Data.Array.Accelerate.IO as AI
import qualified Data.Vector as V 


wignerVille :: A.Array A.DIM1 (ADC.Complex Double) -> V.Vector (A.Array A.DIM1 Double)
wignerVille arr = 
  let vectors = createRVectors arr
      dim1 = (A.Z A.:. (S.size $ shape arr))
  in V.map (ALN.run . (A.map ADC.real) . AMF.fft1D' AMF.Forward dim1) vectors -- `VST.using` (VST.parVector 1)

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

{-# INLINE createRVector #-}
createRVector :: A.Array A.DIM1 (ADC.Complex Double) -> Int -> VS.Vector Int -> VS.Vector Int -> ALN.Acc (A.Array A.DIM1 (ADC.Complex Double))
createRVector arr t inds tau = 
  let sh = shape arr
      leng = {-# SCC leng #-} S.size sh
      shTau = {-# SCC shtau #-}(A.Z A.:. VS.length tau)
      shInds = {-# SCC shinds #-}(A.Z A.:. VS.length inds)
      accArr = {-# SCC accArr #-}A.use arr
      arrInds = {-# SCC arrInds #-}A.use $ AI.fromVectors shInds inds
      expT = {-# SCC expT #-}A.constant (t - 1) 
      accTau = {-# SCC accTau #-}A.use $ AI.fromVectors shTau tau
      accI1 = {-# SCC accI1 #-}A.map ( + expT) accTau
      accI2 = {-# SCC accI2 #-}A.map ((-) (expT)) accTau
      conjAccArr = {-# SCC conjAccArr #-}A.map ADC.conjugate accArr
      v1 = {-# SCC v1 #-}A.scatter accI1 (A.fill (A.constant sh) 0.0) accArr
      v1r = {-# SCC v1r #-}A.reverse $ A.scatter (A.reverse accI1) (A.fill (A.constant sh) 0.0) (A.reverse accArr)
      v2 = {-# SCC v2 #-}A.scatter accI2 (A.fill (A.constant sh) 0.0) conjAccArr
      v2r = {-# SCC v2r #-}A.reverse $ A.scatter (A.reverse accI2) (A.fill (A.constant sh) 0.0) (A.reverse conjAccArr)
      v3 = A.zipWith (*) v1 v2
      v3r = A.zipWith (*) v1r v2r
  in if (t <= leng `mround` 2) 
     then A.scatter arrInds (A.fill (A.constant sh) 0.0) v3
     else A.scatter arrInds (A.fill (A.constant sh) 0.0) v3r
{-
sh :: A.Array A.DIM1 (ADC.Complex Double) -> A.DIM1
sh arr = shape arr

leng :: A.DIM1 -> Int
leng sh = S.size sh

shT :: VS.Vector Int -> A.DIM1
shT tau = (A.Z A.:. VS.length tau)

arrInds :: VS.Vector Int -> ALP.Acc (A.Array A.DIM1 Int)
arrInds inds = 
  let shInds = shT inds
  in A.use $ AI.fromVectors shInds inds

accTau :: VS.Vector Int -> A.Acc (A.Array A.DIM1 Int)
accTau tau = 
  let shTau = shT tau 
  in A.use $ AI.fromVectors shTau tau

expT :: A.Exp Int
expT = A.constant 5

accI1 :: A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int)
accI1 accTau = A.map ( + expT) accTau

accI2 :: A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int)
accI2 accTau =  A.map ((-) expT) accTau

conjAccArr :: A.Acc (A.Array A.DIM1 (ADC.Complex Double)) -> A.Acc (A.Array A.DIM1 (ADC.Complex Double))
conjAccArr accArr = A.map ADC.conjugate accArr

v :: A.Acc (A.Array A.DIM1 Int) -> A.DIM1 -> A.Acc (A.Array A.DIM1 (ADC.Complex Double)) -> A.Acc (A.Array A.DIM1 (ADC.Complex Double))
v accI1 sh accArr = A.scatter accI1 (A.fill (A.constant sh) 0.0) accArr

v3 :: A.Acc (A.Array A.DIM1 (ADC.Complex Double)) -> A.Acc (A.Array A.DIM1 (ADC.Complex Double)) -> A.Acc (A.Array A.DIM1 (ADC.Complex Double))
v3 v1 v2 = A.zipWith (*) v1 v2 -}
{-
v4 :: A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int)
v4 tau arr = 
  let zeros = A.fill (A.shape arr) 0
  in
  A.permute (+) zeros (\ix -> if (ix <= A.shape tau) then A.index1 (tau A.! ix) else ix) arr -}


createParamsStruct :: A.Array A.DIM1 (ADC.Complex Double) -> V.Vector (Int, VS.Vector Int, VS.Vector Int) -- V.Vector (ALN.Acc (A.Array A.DIM1 (ADC.Complex Double)))
createParamsStruct arr = 
  let leng = S.size $ shape arr
      ts = V.enumFromN 1 leng
      taumaxs = V.force $ V.map (taumax leng) ts
      taus = V.force $ V.map tau taumaxs
      indicess = V.force $ V.map (indices leng) taus
  in V.zipWith3 (,,) ts indicess taus

createRVectors :: A.Array A.DIM1 (ADC.Complex Double) -> V.Vector (ALN.Acc (A.Array A.DIM1 (ADC.Complex Double)))
createRVectors arr =V.force $ V.map ( uncurry3 (createRVector arr)) $ createParamsStruct arr

--fftMap :: V.Vector (A.Acc (A.Array A.DIM1 (ADC.Complex Double))) -> V.Vector (A.Acc (A.Array A.DIM1 (ADC.Complex Double)))
-- fftMap 

uncurry3 :: (a -> b -> c -> d) -> ((a, b, c) -> d)
uncurry3 f ~(a,b,c) = f a b c 
{-
ind :: Int -> VS.Vector Int ->  VS.Vector Int
ind leng tau = VS.map (\x -> ((x + leng) `mod` (leng)) + 1) tau

createRVector :: VS.Vector Double -> Int -> R.Array R.D R.DIM2 Double
createRVector xs t = let leng = VS.length xs
                   in R.fromFunction (R.Z R.:. leng R.:. leng) (\(R.Z R.:. x R.:. t) -> create x t leng)
                    where
                    	create x t lengt = createElements (x+1) t lengt xs

createGPU :: R.Array R.D R.DIM2 Double -> A.Acc (A.Array A.DIM0 Double) -> A.Acc (A.Array A.DIM2 (ADC.Complex Double))
createGPU rArr zero = let pp :: (R.Array AI.A R.DIM2 Double)
                          (Just pp) =  (AI.computeAccP rArr)
                      in  A.map (\x -> ADC.mkPolar x (A.the zero)) $ A.use $ AI.fromRepa pp

--createFFT :: A.Acc (A.Array A.DIM1 (ADC.Complex Double)) -> A.Acc (A.Array A.DIM Int) -> A.Acc (A.Array A.DIM1 (ADC.Complex Double))
--createFFT arr leng = let dim1 = (A.Z A.:. leng A.:. leng)
--                     in AMF.fft1D' AMF.Forward dim1 array 

createElements :: Int -> Int -> Int -> VS.Vector Double -> Double
createElements x t leng dat = let !taum = taumax t leng
                                  !tu = tau taum
                                  !indexes = ind leng tu
                                  !half = mround leng 2
                              in case (VS.elemIndex x indexes) of Nothing -> if ((leng - half) >= t)&&(t >(half))&&((x==half+1))
                              	                                             then middleElem half dat t
                              	                                             else 0.0 
                                                                  (Just n) -> createElement (tu VS.! n) dat t

middleElem :: Int -> VS.Vector Double -> Int -> Double
middleElem half xs t = (xs VS.! (t + half - 1))*(xs VS.! (t - half -1))


createElement :: Int -> VS.Vector Double -> Int -> Double
createElement tau xs time = (xs VS.! (time + tau))*(xs VS.! (time - tau))
-}
mround :: Int -> Int -> Int 
mround a !b = let half = b `div` 2
             in if (a `mod` b) >= half
                then ((a `div` b) + 1)
                else (a `div` b)
