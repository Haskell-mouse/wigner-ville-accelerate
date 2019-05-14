{-# LANGUAGE FlexibleContexts, ScopedTypeVariables #-}
-- |
-- Module      : Data.Array.Accelerate.Math.Choi-Williams
-- Copyright   : [2017] Rinat Stryungis
-- License     : BSD3
--
-- Maintainer  : Rinat Stryungis <lazybonesxp@gmail.com>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
-- Functions, needed for computation of Choi-Williams transform using the accelerate-fft library.
-- Original alghorithm : https://calhoun.nps.edu/bitstream/handle/10945/4422/09Dec_Hollinger.pdf
--
-- This module uses the accelerate-fft library. And the base implementation of fft
-- uses a naive divide-and-conquer fft implementation
-- whose absolute performance is appalling. It also requires that you know on
-- the Haskell side the size of the data being transformed, and that this is
-- a power-of-two in each dimension.
--
-- For performance, compile accelerate-fft against the foreign library bindings (using any
-- number of '-fllvm-ptx', and '-fllvm-cpu' for the accelerate-llvm-ptx, and
-- accelerate-llvm-native backends, respectively), which have none of the above
-- restrictions.
-- Both of this flags are enabled by default.

module Data.Array.Accelerate.Math.BornJordan where

import Data.Array.Accelerate as A
import Data.Array.Accelerate.Control.Lens
import qualified Data.Array.Accelerate.Data.Complex as ADC
import qualified Data.Array.Accelerate.Math.FFT as AMF
import Data.Array.Accelerate.Math.Hilbert
import qualified Data.Array.Accelerate.Math.PseudoWigner as P
import Data.Array.Accelerate.Math.Wigner'

-- | Create matrix of x (mu-tau)*(*(mu+tau)). And replicate it to third dimension.

createMatrixHalf1 :: (Elt (ADC.Complex e),A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => Acc (Array DIM1 (ADC.Complex e))   -- ^ Data array
  -> Acc (Array DIM1 Int)               -- ^ taumaxs
  -> Acc (Array DIM1 Int)               -- ^ limits
  -> Exp Int
  -> Exp Int
  -> Acc (Array DIM3 (ADC.Complex e))
createMatrixHalf1 arr taumaxs lims uWindow nWindow = matrixQ  -- | Add the time dimension
  where
    matrixF = createMatrixHalf arr taumaxs lims
    matrixQ = A.generate (A.index3 uLeng nLeng leng)
                (\sh -> let (A.Z A.:. mu A.:. tau A.:. t) = A.unlift sh
                            muMin = A.min 0 (A.abs (t - uLengH))
                            muMax = A.min (leng - 1) (t + uLengH)
                            newMu = mu + (t - uLengH )
                        in A.cond ((newMu A.>= muMin A.|| (t A.> uLengH)) A.&& newMu A.<= muMax)
                                  (matrixF A.! (A.index2 newMu tau))
                                  (A.constant $ 0.0 ADC.:+ 0.0))

    leng = A.length arr
    uLeng = A.min uWindow leng
    uLengH = uLeng `A.div` 2
    nLeng = ((nWindow `A.div` 2) + 1)
    movedUpTau tau taum = A.abs $ (taum + tau) `A.mod` leng
    gen mu tau lim taum = A.cond ((movedUpTau tau taum) A.< lim) (generateValue arr mu ((movedUpTau tau taum) - taum)) 0

-- | Create 3D array with kernel function ()

coreFunctionHalf :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)  --
  => Exp Int                            -- ^ length in each dimention
  -> Exp e                              -- ^ sigma
  -> Acc (Array DIM1 Int)               -- ^ taumaxs
  -> Acc (Array DIM1 Int)               -- ^ limits
  -> Maybe (Acc (Array DIM2 e))
  -> Exp Int                            --
  -> Exp Int                            --
  -> A.Acc (Array DIM2 e)
coreFunctionHalf leng sigma taumaxs lims mWindowArray uWindow nWindow =
  let uLeng = A.min uWindow leng
      nLeng = (nWindow `A.div` 2) + 1
      applyWindowArray =
        case mWindowArray of
        Nothing -> id
        Just windowArray -> A.zipWith (*) windowArray
      unnormalised = A.generate (A.index2 uLeng nLeng) (\sh ->
                        let (A.Z A.:.mu A.:. tau) = A.unlift sh
                            lim = lims A.!! newMu
                            taum = taumaxs A.!! newMu
                            newMu = mu - (uLeng `A.div` 2)
                        in genCore newMu tau sigma nWindow)
      in applyWindowArray unnormalised

coreFunctionHalf3D :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)  --
  => Exp Int                            -- ^ length in each dimention
  -> Exp e                              -- ^ sigma
  -> Acc (Array DIM1 Int)               -- ^ taumaxs
  -> Acc (Array DIM1 Int)               -- ^ limits
  -> Maybe (Acc (Array DIM2 e))
  -> Exp Int                            --
  -> Exp Int                            --
  -> A.Acc (Array DIM3 e)
coreFunctionHalf3D leng sigma taumaxs lims mWindowArray uWindow nWindow =
  let core = coreFunctionHalf leng sigma taumaxs lims mWindowArray uWindow nWindow
      nLeng = ((nWindow `A.div` 2) + 1)
      uLeng = A.min uWindow leng
      uLengH = uLeng `A.div` 2
  in  A.generate (A.index3 uLeng nLeng leng)
                 (\sh -> let (A.Z A.:. mu A.:. tau A.:. t) = A.unlift sh
                         in A.cond ((t A.< uLengH A.&& mu A.< uLengH - t) A.|| (t A.>= leng - uLengH A.&& mu A.>= (uLengH + (leng - t)))) --((newMu A.>= muMin A.|| (t A.> uLengH)) A.&& newMu A.<= muMax)
                                    (A.constant 0.0)
                                    (core A.! (A.index2 mu tau)))

-- | generate each value of 3D array with kernel.

genCore :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e, A.Eq e, A.FromIntegral Int e)
  => Exp Int                            -- ^ mu
  -> Exp Int
  -> Exp e                              -- ^ sigma
  -> Exp Int
  -> Exp e
genCore mu tau sigma nWindow =
  let mu1 = (A.fromIntegral mu) -- -uWindow/2 < u < uWindow/2 - 1
      tau1 = A.fromIntegral tau
  in A.cond (A.abs mu1 A.> A.abs tau) 0.0 (1.0/((A.abs (tau1*2)) + 1.0 ))

-- | Product element-by-element result of amatrix and result of coreFunction, and fold it over mu.

sFunc :: (Elt (ADC.Complex e),AMF.Numeric e,A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => A.Acc (Array DIM3 e)               -- ^ core 3D array
  -> A.Acc (Array DIM3 (ADC.Complex e)) -- ^ result of amatrix
  -> A.Acc (Array DIM2 (ADC.Complex e))
sFunc core aMatrix =  A.fold (+) 0 $ A.zipWith (*) (A.map makeComplex $ transposeOn _1 _3 core) (transposeOn _1 _3 aMatrix)


summedOverMu :: (Elt (ADC.Complex e), A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e, AMF.Numeric e, Prelude.Fractional (A.Exp e))
  => Acc (Array DIM1 (ADC.Complex e))   -- ^ Data array
  -> Acc (Array DIM1 Int)               -- ^ taumaxs
  -> Acc (Array DIM1 Int)               -- ^ limits
  -> Exp e                              -- ^ sigma
  -> Maybe (Acc (Array DIM1 e), Acc (Array DIM1 e))
  -> Exp Int                 -- ^ uWindow
  -> Exp Int                 -- ^ nWindow
  -> Acc (Array DIM2 (ADC.Complex e))
summedOverMu arr taumaxs lims sigma  mWindowArrays uWindow nWindow =
  let leng = A.length arr
      nLeng = nWindowHalf + 1
      hlf = (leng `A.div` 2) - 1
      nWindowHalf = nWindow `A.div` 2

      mWindowArray = case mWindowArrays of
        Nothing -> Nothing
        Just (uWindowArray,nWindowArray) -> Just $ A.zipWith (*) uWindowArray2D nWindowArray2D
                                            where  uWindowArray2D = A.replicate (A.lift $ A.Z A.:. All A.:. nLeng) uWindowArray
                                                   nWindowArrayHalf = A.drop nWindowHalf nWindowArray
                                                   nWindowArray2D = A.transpose $ A.replicate (A.lift $ A.Z A.:. All A.:. uWindow) nWindowArrayHalf

      zeroesForAmatrix = A.fill (A.index2 leng ((leng `A.div` 2) - nLeng)) (A.constant $ 0.0 ADC.:+ 0.0)
      amatrix = (createMatrixHalf1 arr taumaxs lims uWindow nWindow)
      coreMatrix = coreFunctionHalf leng sigma taumaxs lims mWindowArray uWindow nWindow
      coreMatrixNormalised = normalise uWindow nLeng leng $
        coreFunctionHalf3D leng sigma taumaxs lims mWindowArray uWindow nWindow
      coreMatrix3D = coreMatrixNormalised

      summed = (sFunc coreMatrix3D amatrix) A.++ zeroesForAmatrix
  in generate (A.index2 leng leng) (\sh -> let (A.Z A.:.time A.:. tau) = A.unlift sh
                                           in cond (tau A.< hlf) (summed A.! (A.index2 time tau)) (ADC.conjugate $ summed A.! (A.index2 time (leng - tau))))


normalise :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => A.Exp Int
  -> A.Exp Int
  -> A.Exp Int
  -> A.Acc (Array DIM3 e)
  -> A.Acc (Array DIM3 e)
normalise lengMu lengTau leng date =
  let transposed = transposeOn _1 _3 date
      summedColumn y x = A.sfoldl (+) 0.0 (A.lift (A.Z A.:. (y :: A.Exp Int) A.:. (x :: A.Exp Int) )) transposed
      summed = A.generate (A.index2 leng lengTau) (\sh ->
                  let (A.Z A.:. (t :: A.Exp Int) A.:. (tau :: A.Exp Int)) = A.unlift sh
                  in summedColumn t tau)
      matrix3 = A.replicate (A.lift $ A.Z A.:.All A.:. A.All A.:.lengMu) summed
  in transposeOn _1 _3 $ A.zipWith (/) transposed matrix3

{-

genCore :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => Exp Int                            -- ^ n - lag
  -> Exp Int                            -- ^ l - time
  -> Exp Int                            -- ^ mu
  -> Exp Int                            -- ^ length of array
  -> Exp Int
  -> Exp Int
  -> Exp e
genCore mu tau time leng uWindow nWindow =
  let mu1 = ((A.fromIntegral mu) - udiff )  :: A.Exp Int          -- -N/2 < u < N/2 - 1
      time1 = (A.fromIntegral time) - h :: A.Exp Int             -- The same limits
      tau1 = A.fromIntegral tau -- cond ((A.fromIntegral n) A.< h) (A.fromIntegral n) (A.fromIntegral (leng - n)) :: A.Exp Int
      h = A.fromIntegral (leng `div` 2) :: A.Exp Int
      udiff = A.min (A.fromIntegral $ (uWindow `div` 2) + 1) h
  in (cond
       (mu1 A.<= (time1) + (A.abs tau1) A.&&
        mu1 A.>= (time1) - (A.abs tau1))
          (cond (tau1 A.== 0 A.&& mu1 A.== 0) 1.0 (1.0/((A.fromIntegral $ (A.abs (tau1*2))) + 1.0 ))) 0.0)

-}
