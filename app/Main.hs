{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE OverloadedStrings #-}

module Main where

import RwData
import Wigner
import Hilbert
import Data.Either
import qualified Data.Array.IArray as IA
import qualified Data.Array.Accelerate as A
import qualified Data.Array.Accelerate.Interpreter as ALI
import qualified Data.Array.Accelerate.LLVM.Native as ALN
import qualified Data.Array.Accelerate.LLVM.PTX as ALP
import Data.Array.Accelerate.Array.Sugar as S 
import System.Environment
import Data.Monoid
import Debug.Trace 
import System.IO
import GHC.Float
import TextShow.Data.Array
import qualified TextShow as TS 
import qualified Data.Array.Accelerate.IO as AI
import qualified Data.Array.Accelerate.Data.Complex as ADC
import qualified Data.Vector.Storable as VS
import qualified Data.Text.IO as TI
import qualified Data.Text.Lazy.IO as TLI
import qualified Data.Double.Conversion.Text as DT 
import qualified Data.Text as T
import qualified Control.Parallel.Strategies as CP
import Data.List.Split
import qualified Data.Vector as V


main :: IO ()
main = do
  args <- getArgs
  let num = A.fromList (A.Z) $ [w]
      fr = S.toList $ ALP.run1 makeGlobal num
      w = read $ head args :: Int 
  file <- openFile "wigner.txt" WriteMode
  writeString file fr w 1
  hClose file
  
writeString :: Handle -> [Float] -> Int -> Int -> IO ()
writeString _ [] _ _  = return ()
writeString file (x:xs) w n = do 
  TI.hPutStr file $ DT.toPrecision 5 (float2Double x) <> " "
  if (n `mod` w) == 0 
  then do 
    TI.hPutStr file "\n"
    writeString file xs w (n + 1)
  else writeString file xs w (n + 1)

readLines :: FilePath -> IO [String]
readLines = fmap lines . readFile

makeDouble :: [String] -> [Double]
makeDouble = map read

--start :: [A.Array A.DIM1 Double] -> Period -> Maybe (A.Array A.DIM2 (ADC.Complex Double))
--start vectors period = Just  (wignerVille (ALI.run $ hylbert vectors))
