{-# LANGUAGE BangPatterns #-}

module Main where

import RwData
import Wigner
import Hilbert
import Data.Either
import qualified Data.Array.Accelerate as A
import qualified Data.Array.Accelerate.Interpreter as ALI
import qualified Data.Array.Accelerate.LLVM.Native as ALN
--import qualified Data.Array.Accelerate.LLVM.PTX as ALP
import Data.Array.Accelerate.Array.Sugar as S 
import System.Environment
import qualified Data.Array.Accelerate.IO as AI
import qualified Data.Array.Accelerate.Data.Complex as ADC
import qualified Data.Vector.Storable as VS
import qualified Data.Text.IO as TI
import qualified Data.Double.Conversion.Text as DT 
import qualified Data.Text as T 
import qualified Data.Vector as V


main :: IO ()
main = do
  args <- getArgs
  let fname = head args
  strNum <- readLines fname
  let nums = makeDouble strNum
  let leng = length nums
      vect = A.fromList (A.Z A.:. leng) nums
      v = ALN.run $ hylbert vect
      z = wignerVille v
      t = V.map AI.toVectors z
      vToText = VS.foldr (\x buf -> T.append buf (T.append (DT.toShortest x) (T.pack " "))) T.empty
      vsToText = V.foldr (\x buf -> T.append buf x) T.empty . V.map ((flip T.append (T.pack "\n")) . vToText)
      text = vsToText t
  TI.writeFile "wigner.txt" text


readLines :: FilePath -> IO [String]
readLines = fmap lines . readFile

makeDouble :: [String] -> [Double]
makeDouble = map read

--start :: [A.Array A.DIM1 Double] -> Period -> Maybe (A.Array A.DIM2 (ADC.Complex Double))
--start vectors period = Just  (wignerVille (ALI.run $ hylbert vectors))
