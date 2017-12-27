{-# LANGUAGE FlexibleContexts #-}
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
-- Original alghorithm : http://www.dtic.mil/dtic/tr/fulltext/u2/a514435.pdf
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

import Data.Array.Accelerate.Math.Wigner'
import Data.Array.Accelerate.Math.Hilbert
import qualified Data.Array.Accelerate.Math.PseudoWigner as P
import Data.Array.Accelerate as A
import qualified Data.Array.Accelerate.Math.FFT as AMF
import qualified Data.Array.Accelerate.Data.Complex as ADC


-- | Create matrix of x (t-tau)*(*(x+tau)). And replicate it to third dimension. 

amatrix :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e) 
  => Acc (Array DIM1 (ADC.Complex e))   -- ^ Data array
  -> Acc (Array DIM1 Int)               -- ^ taumaxs
  -> Acc (Array DIM1 Int)               -- ^ limits 
  -> Acc (Array DIM3 (ADC.Complex e))
amatrix arr taumx lims = 
  let a = A.transpose $ createMatrix arr taumx lims
      leng = A.length arr  
  in A.replicate (A.lift $ A.Z A.:. All A.:. leng A.:. All) a

amatrix2 :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e) 
  => Acc (Array DIM1 (ADC.Complex e))   -- ^ Data array
  -> Acc (Array DIM1 Int)               -- ^ taumaxs
  -> Acc (Array DIM1 Int)               -- ^ limits 
  -> Acc (Array DIM2 (ADC.Complex e))
amatrix2 arr taumx lims = createMatrix arr taumx lims

amatrix_w :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e) 
  => Acc (Array DIM1 (ADC.Complex e))   -- ^ Data array
  -> Acc (Array DIM1 Int)               -- ^ taumaxs
  -> Acc (Array DIM1 Int)               -- ^ limits 
  -> Acc (Array DIM1 e)                 -- ^ Smoothing window. Length of it must be odd.
  -> Acc (Array DIM3 (ADC.Complex e))
amatrix_w arr taumx lims window = 
  let a = A.transpose $ P.createMatrix arr window taumx lims
      leng = A.length arr  
  in A.replicate (A.lift $ A.Z A.:. All A.:. leng A.:. All) a

-- | Create 3D array with kernel function ()

coreFunction :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)  -- 
  => Exp Int                            -- ^ length in each dimention
  -> Exp e                              -- ^ alpha 
  -> A.Acc (Array DIM3 e)
coreFunction leng alpha =  
  A.generate (A.index3 leng leng leng) (\sh -> 
                let (A.Z A.:.u A.:. l A.:. n) = A.unlift sh 
                in genCore u l n leng alpha)

-- | generate each value of 3D array with kernel. 

genCore :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e) 
  => Exp Int                            -- ^ n 
  -> Exp Int                            -- ^ l 
  -> Exp Int                            -- ^ mu 
  -> Exp Int                            -- ^ length of array
  -> Exp e                              -- ^ alpha
  -> Exp e
genCore n l u leng alpha = 
  let u1 = ((A.fromIntegral u) - h )           -- -N/2 < u < N/2 - 1
      l1 = (A.fromIntegral l) - h              -- The same limits
      n1 = cond ((A.fromIntegral n) A.< h) (A.fromIntegral n) (A.fromIntegral (n - leng)) 
      h = A.fromIntegral (leng `div` 2)
  in cond (u1 A.<= l1 + alpha*(A.abs n1) A.&& u1 A.>= l1 - alpha*(A.abs n1)) (cond (n1 A.== 0) (cond (u1 A.== l1) 1.0 0.0) (1.0/(2.0*alpha*(A.abs n1)))) 0.0

-- | Product element-by-element result of amatrix and result of coreFunction, and fold it over mu. 

sFunc :: (A.RealFloat e, A.IsFloating e, A.FromIntegral Int e, Elt e)
  => A.Acc (Array DIM3 e)               -- ^ core 3D array       
  -> A.Acc (Array DIM3 (ADC.Complex e)) -- ^ result of amatrix
  -> A.Acc (Array DIM2 (ADC.Complex e))
sFunc core aMatrix =  A.fold (+) 0 $ A.zipWith (*) (A.map makeComplex core) aMatrix