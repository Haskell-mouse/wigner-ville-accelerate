{-# LANGUAGE FlexibleContexts, TypeFamilies #-}

-- |
-- Module      : Data.Array.Accelerate.Math.Wigner`
-- Copyright   : [2017] Rinat Stryungis
-- License     : BSD3
--
-- Maintainer  : Rinat Stryungis <lazybonesxp@gmail.com>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
-- Helper functions for a Wigner transform using the accelerate-fft library.

module Data.Array.Accelerate.Math.Wigner' where

import qualified Data.Array.Accelerate as A
import Data.Array.Accelerate.Array.Sugar as S
import qualified Data.Array.Accelerate.Data.Complex as ADC
import qualified Data.Array.Accelerate.Math.FFT as AMF

taumax :: A.Exp Int -> A.Exp Int -> A.Exp Int
taumax leng t = min (min t (leng - t - 1) ) (A.round (((A.fromIntegral leng)/2.0) - 1.0 :: A.Exp Double))

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

moveUp :: A.Acc (A.Array A.DIM1 Int) -> A.Exp Int -> A.Exp DIM2 -> A.Exp DIM2
moveUp taumaxs leng sh =
  let taum t = taumaxs A.!! t
  in (\(t, tau) -> A.index2 t ((tau + (taum t)) `A.mod` leng)) $ A.unlift $ A.unindex2 sh

generateValue :: (A.RealFloat e, Elt (ADC.Complex e), Fractional (A.Exp e), Floating (A.Exp e), A.IsFloating e, A.FromIntegral Int e, Elt e) => A.Acc (A.Array A.DIM1 (ADC.Complex e)) -> A.Exp Int -> A.Exp Int -> A.Exp (ADC.Complex e)
generateValue arr time tau = (arr A.!! (time + tau)) * (ADC.conjugate $ arr A.!! (time - tau))


createMatrix :: (Elt (ADC.Complex e), A.RealFloat e, Fractional (A.Exp e), Floating (A.Exp e), A.IsFloating e, A.FromIntegral Int e, Elt e) => A.Acc (A.Array A.DIM1 (ADC.Complex e)) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM2 (ADC.Complex e))
createMatrix arr taumaxs lims = raw
  where
    raw = A.generate (A.index2 leng leng) (\sh -> let (A.Z A.:. t A.:. tau) = A.unlift sh
                                                      lim = lims A.!! t
                                                      taum = taumaxs A.!! t
                                                  in gen t tau lim taum)
    leng = A.length arr
    movedUpTau tau taum = A.abs $ (taum + tau) `A.mod` leng
    gen t tau lim taum = A.cond ((movedUpTau tau taum) A.< lim) (generateValue arr t ((movedUpTau tau taum) - taum)) 0

createMatrixHalf :: (Elt (ADC.Complex e), A.RealFloat e, Fractional (A.Exp e), Floating (A.Exp e), A.IsFloating e, A.FromIntegral Int e, Elt e) => A.Acc (A.Array A.DIM1 (ADC.Complex e)) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM2 (ADC.Complex e))
createMatrixHalf arr taumaxs lims = {-A.backpermute (A.index2 leng leng) (moveUp taumaxs leng) $ -}raw
  where
    raw = A.generate (A.index2 leng ((leng `A.div` 2) + 1)) (\sh -> let (A.Z A.:. t A.:. tau) = A.unlift sh
                                                                        lim = lims A.!! t
                                                                        taum = taumaxs A.!! t
                                                                    in gen t tau lim taum)
    leng = A.length arr
    movedUpTau tau taum = A.abs $ (taum + tau) `A.mod` leng
    gen t tau lim taum = A.cond ((movedUpTau tau taum) A.< lim) (generateValue arr t ((movedUpTau tau taum) - taum)) 0

createMatrixFull :: (Elt (ADC.Complex e), A.RealFloat e, Fractional (A.Exp e), Floating (A.Exp e), A.IsFloating e, A.FromIntegral Int e, Elt e) => A.Acc (A.Array A.DIM1 (ADC.Complex e)) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM1 Int) -> A.Acc (A.Array A.DIM2 (ADC.Complex e))
createMatrixFull arr taumaxs lims =
  let hlfArr = createMatrixHalf arr taumaxs lims
      leng = A.length arr
      hlf = (leng `A.div` 2) + 1
  in A.generate (A.index2 leng leng) (\sh -> let (A.Z A.:. t A.:. tau) = A.unlift sh
                                             in A.cond (tau A.< hlf) (hlfArr A.! (A.index2 t tau)) (ADC.conjugate $ (hlfArr A.! (A.index2 t (leng - tau)))))
