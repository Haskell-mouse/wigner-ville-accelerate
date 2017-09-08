{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts#-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE EmptyDataDecls #-}

module Smoothed(WindowFunc(..), spWignerVille, matProduct, mkSmMatrix) where

import Hilbert
import PseudoWigner
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

spWignerVille :: (A.RealFloat e, A.IsFloating e, Elt e, A.FromIntegral Int e) => 
  ((A.Acc (A.Array A.DIM1 e), A.Acc (A.Array A.DIM1 e)), A.Acc (A.Array A.DIM1 (ADC.Complex e))) -> A.Acc (A.Array A.DIM2 e)
spWignerVille ((fWindow, tWindow), arr) = 
  let pwv = pWignerVille (fWindow,arr)
      leng = A.unit $ A.length arr
      tleng = A.length tWindow
      smoothMatrix = mkSmMatrix leng tWindow
  in A.map (\x -> x/(A.fromIntegral tleng)) $ matProduct pwv smoothMatrix

mkSmMatrix :: (A.RealFloat e, Elt e) => 
  A.Acc (A.Scalar Int) -> A.Acc (A.Array A.DIM1 e) -> A.Acc (A.Array A.DIM2 e)
mkSmMatrix arrLeng window = 
  A.generate (A.index2 leng leng) (\sh -> let (A.Z A.:. x A.:. t) = A.unlift sh in gen x t)
  where 
    leng = A.the arrLeng
    wleng = A.length window
    gen x t = A.cond ((x A.< t - (wleng - 1) `A.div` 2) A.|| (x A.> t + (wleng - 1) `A.div` 2)) 0 (window A.!! (x + (wleng `A.div` 2) - t))

matProduct :: (A.RealFloat e, Elt e) => 
  A.Acc (A.Array A.DIM2 e) -> A.Acc (A.Array A.DIM2 e) -> A.Acc (A.Array A.DIM2 e)
matProduct wign smooth = 
  let A.Z A.:. wx A.:. _ = A.unlift $ A.shape wign :: A.Z A.:. A.Exp Int :. A.Exp Int 
      A.Z A.:. _ A.:. sy = A.unlift $ A.shape smooth :: A.Z A.:. A.Exp Int :. A.Exp Int
      aRep = A.replicate (A.lift $ (A.Z A.:. A.All A.:. sy A.:. A.All)) wign
      bRep = A.replicate (A.lift $ (A.Z A.:. wx A.:. A.All A.:. A.All)) (A.transpose smooth)
  in A.fold (+) 0 $ A.zipWith (*) aRep bRep 