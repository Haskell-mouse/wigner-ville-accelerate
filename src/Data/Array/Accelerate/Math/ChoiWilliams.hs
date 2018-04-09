{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables#-}
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

module Data.Array.Accelerate.Math.ChoiWilliams where

import Data.Array.Accelerate.Math.Wigner'
import Data.Array.Accelerate.Math.Hilbert
import qualified Data.Array.Accelerate.Math.PseudoWigner as P
import Data.Array.Accelerate as A
import qualified Data.Array.Accelerate.Math.FFT as AMF
import qualified Data.Array.Accelerate.Data.Complex as ADC


-- | Create matrix of x (t-tau)*(*(x+tau)). And replicate it to third dimension. 

amatrix :: (Elt (ADC.Complex e),A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e) 
  => Acc (Array DIM1 (ADC.Complex e))   -- ^ Data array
  -> Acc (Array DIM1 Int)               -- ^ taumaxs
  -> Acc (Array DIM1 Int)               -- ^ limits 
  -> Acc (Array DIM3 (ADC.Complex e))
amatrix arr taumx lims = 
  let a = A.transpose $ createMatrix arr taumx lims
      leng = A.length arr  
  in A.replicate (A.lift $ A.Z A.:. All A.:. leng A.:. All) a

amatrix_w :: (A.RealFloat e,Elt (ADC.Complex e),AMF.Numeric e, A.IsFloating e, A.FromIntegral Int e, Elt e) 
  => Acc (Array DIM1 (ADC.Complex e))   -- ^ Data array
  -> Acc (Array DIM1 Int)               -- ^ taumaxs
  -> Acc (Array DIM1 Int)               -- ^ limits 
  -> Acc (Array DIM1 e)                 -- ^ Smoothing window. Length of it must be odd.
  -> Acc (Array DIM3 (ADC.Complex e))
amatrix_w arr taumx lims window = 
  let a = A.transpose $ P.createMatrix arr window taumx lims
      leng = A.length arr  
  in A.replicate (A.lift $ A.Z A.:. All A.:. All A.:. leng) a

-- | Create 3D array with kernel function ()

coreFunction :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)  -- 
  => Exp Int                            -- ^ length in each dimention
  -> Exp e                              -- ^ sigma 
  -> A.Acc (Array DIM3 e)
coreFunction leng sigma =  
  A.generate (A.index3 leng leng leng) (\sh -> 
                let (A.Z A.:.u A.:. l A.:. n) = A.unlift sh 
                in genCore u l n leng sigma)

coreFunction2 :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)  -- 
  => Exp Int                            -- ^ length in each dimention
  -> Exp e                              -- ^ sigma 
  -> A.Acc (Array DIM3 e)
coreFunction2 leng sigma =  
 (flip testfunc) leng $ A.generate (A.index3 leng leng leng) (\sh -> 
                            let (A.Z A.:.u A.:. l A.:. n) = A.unlift sh 
                            in genCore u l n leng sigma)   -- reverse n l and u

-- | generate each value of 3D array with kernel. 

genCore :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e) 
  => Exp Int                            -- ^ n - Time lag 
  -> Exp Int                            -- ^ l - Time
  -> Exp Int                            -- ^ mu - 
  -> Exp Int                            -- ^ length of array
  -> Exp e                              -- ^ sigma
  -> Exp e
genCore n l u leng sigma = 
  let u1 = ((A.fromIntegral u) - h )           -- -N/2 < u < N/2 - 1
      l1 = (A.fromIntegral l) - h              -- The same limits
      n1 = cond ((A.fromIntegral n) A.< h) (A.fromIntegral n) (A.fromIntegral (n - leng)) 
      h = A.fromIntegral (leng `div` 2)
  in cond (n1 A.== 0) (cond (u1 A.== l1) 1.0 0.0) (((sqrt(sigma))/sqrt ((4.0*pi*n1*n1))) * exp ((-1.0)*(((u1 - l1)*(u1 - l1))*sigma/(4.0*n1*n1)) ) )

-- | Product element-by-element result of amatrix and result of coreFunction, and fold it over mu. 

sFunc :: (Elt (ADC.Complex e),AMF.Numeric e,A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => A.Acc (Array DIM3 e)               -- ^ core 3D array       
  -> A.Acc (Array DIM3 (ADC.Complex e)) -- ^ result of amatrix
  -> A.Acc (Array DIM2 (ADC.Complex e))
sFunc core aMatrix =  A.fold (+) 0 $ A.zipWith (*) (A.map makeComplex core) aMatrix

testfunc :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => A.Acc (Array DIM3 e)
  -> A.Exp Int
  -> A.Acc (Array DIM3 e)
testfunc date leng = 
  let gentest x y = A.sfoldl (+) 0.0 (A.lift (A.Z A.:. (x :: A.Exp Int) A.:. (y :: A.Exp Int))) date
      matrix2 = A.generate (A.index2 leng leng) (\sh -> 
                let (A.Z A.:.(u :: A.Exp Int) A.:. (l :: A.Exp Int)) = A.unlift sh 
                in gentest u l)
      matrix3 = A.replicate (A.lift $ A.Z A.:.All A.:.All A.:.leng) matrix2
  in A.zipWith (/) date matrix3
