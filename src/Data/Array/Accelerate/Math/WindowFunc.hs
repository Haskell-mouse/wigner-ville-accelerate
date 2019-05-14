{-# LANGUAGE DeriveDataTypeable, FlexibleContexts #-}

-- |
-- Module      : Data.Array.Accelerate.Math.WindowFunc
-- Copyright   : [2017] Rinat Stryungis
-- License     : BSD3
--
-- Maintainer  : Rinat Stryungis <lazybonesxp@gmail.com>
-- Stability   : experimental
-- Portability : non-portable (GHC extensions)
--
-- Creation of window for smoothing in frequency-domain in Pseudo-Wigner-Ville distribuition

module Data.Array.Accelerate.Math.WindowFunc(WindowFunc(..),makeWindow) where

import qualified Data.Array.Accelerate as A
import Data.Data
import Data.Typeable

-- | Function of the window. Rect - Rectangle.

data WindowFunc = Rect | Sin | Lanczos | Hanning | Hamming | Bartlett
  deriving (Read, Show, Data, Typeable, Eq)

-- | Creates new window (1D array of odd length) with length and window function.
-- For example
--   win1 = makeWindow Sin lentgh
-- Where length has type Acc (Scalar Int)

makeWindow :: (A.RealFloat e, Fractional (A.Exp e), Floating (A.Exp e), A.IsFloating e, A.FromIntegral Int e, Ord e) =>
  WindowFunc -> A.Acc (A.Scalar Int) -> A.Acc (A.Array A.DIM1 e)
makeWindow func leng =
  let gen = A.generate (A.index1 $ A.the leng)
  in case func of
       Rect -> A.fill (A.index1 $ A.the leng) 1.0
       Sin  -> gen (\sh -> let (A.Z A.:.x) = A.unlift sh in sin (pi*(A.fromIntegral x)/(A.fromIntegral $ A.the leng - 1)))
       Lanczos -> gen (\sh -> let (A.Z A.:.x) = A.unlift sh in sinc ((2*(A.fromIntegral x)/(A.fromIntegral $ A.the leng - 1)) - 1.0))
       Hanning -> gen (\sh -> let (A.Z A.:.x) = A.unlift sh in 0.5 - (0.5 * (cos (2*pi*(A.fromIntegral (x + 1))/(A.fromIntegral $ A.the leng + 1)))))
       Hamming -> gen (\sh -> let (A.Z A.:.x) = A.unlift sh in 0.54 - (0.46 * (cos (2*pi*(A.fromIntegral (x + 1))/(A.fromIntegral $ A.the leng + 1)))))
       Bartlett -> gen (\sh -> let (A.Z A.:.x) = A.unlift sh in 1.0 - A.abs (((A.fromIntegral x)/(A.fromIntegral (A.the leng - 1)/2.0)) - 1.0))

sinc :: (Floating (A.Exp e), A.Elt e, A.Ord e) => A.Exp e -> A.Exp e
sinc x =
  A.cond (ax A.< eps_0) 1 (A.cond (ax A.< eps_2) (1 - x2/6) (A.cond (ax A.< eps_4) (1 - x2/6 + x2*x2/120) ((A.sin x)/x)))
  where
    ax = A.abs x
    x2 = x*x
    eps_0 = 1.8250120749944284e-8 -- sqrt (6ε/4)
    eps_2 = 1.4284346431400855e-4 --   (30ε)**(1/4) / 2
    eps_4 = 4.043633626430947e-3  -- (1206ε)**(1/6) / 2
