{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE OverloadedStrings #-}

module Main where

import ParseArgs
import ParseFile
import Wigner
import Hilbert
import PseudoWigner
--import Smoothed
import Data.Either
import Data.Attoparsec.Text
import System.Directory
import qualified Data.Array.IArray as IA
import qualified Data.Array.Accelerate as A
import qualified Data.Array.Accelerate.Interpreter as ALI
import qualified Data.Array.Accelerate.LLVM.Native as ALN
import qualified Data.Array.Accelerate.LLVM.PTX as ALP
import Data.Array.Accelerate.Array.Sugar as S 
import System.Environment
import qualified Data.Vector.Split as DS
import Data.Monoid
import Debug.Trace 
import Control.DeepSeq
import System.IO
import GHC.Float
import Data.List
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
  args <- opts
  case args of 
       (OptsWV _ _ _) -> makeWVall args
       (OptsPWV _ _ _ _ _) -> makePWVall args
  
writeString :: Handle -> [Double] -> Int -> Int -> IO ()
writeString _ [] _ _  = return ()
writeString file (x:xs) w n = do 
  TI.hPutStr file $ DT.toPrecision 5 x <> " "
  if (n `mod` w) == 0 
  then do 
    TI.hPutStr file "\n"
    writeString file xs w (n + 1)
  else writeString file xs w (n + 1)  

{-writeString :: Handle -> [Double] -> Int -> Int -> IO ()
writeString file xs w n = do 
  let txt = CP.withStrategy (CP.parListChunk 4000 CP.rdeepseq) $  map (\y -> DT.toPrecision 5 y <> " ") xs
  writeToFile file txt w n
  
{-# INLINE writeToFile #-}
writeToFile :: Handle -> [T.Text] -> Int -> Int -> IO ()
writeToFile _ [] _ _  = return ()
writeToFile file (x:xs) w n = do 
  TI.hPutStrLn file x
  if (n `mod` w) == 0 
  then do 
    TI.hPutStr file "\n"
    writeToFile file xs w (n + 1)
  else writeToFile file xs w (n + 1) -}


makePWVall :: Opts -> IO ()
makePWVall (OptsPWV path dev wlen wfun nosAVGflag) = do
  if (even wlen)
  then error "Window length maust be odd !!"
  else do  
    files <- listDirectory path
    setCurrentDirectory path 
    let fFiles = filter (isSuffixOf ".txt") files
    let window = makeWindow wfun (A.unit $ A.constant wlen) :: A.Acc (A.Array A.DIM1 Double)
        sAVGflag = A.constant $ not nosAVGflag
    mapM_ (makePWV dev window sAVGflag) fFiles

makePWV :: CalcDev -> A.Acc (A.Array A.DIM1 Double) -> A.Exp Bool -> FilePath -> IO ()
makePWV dev window sAVGflag file = do 
  putStrLn $ "processing " ++ file ++ "..."
  text <- TI.readFile file 
  let parseRes = parseOnly parseFile text 
  case parseRes of 
    (Left errStr) -> error $ errStr 
    (Right dataF) -> mapM_ (startPWV dev window file sAVGflag) dataF 

startPWV :: CalcDev -> A.Acc (A.Array A.DIM1 Double) -> FilePath -> A.Exp Bool -> (T.Text,[Double]) -> IO ()
startPWV dev window oldName sAVGflag (file,dataF) = do
  let newFName = oldName ++ "-" ++ T.unpack file ++ ".txt"
  putStr $ "   Creating file " ++ newFName ++ " ..."
  let leng = length dataF
      pData = A.fromList (A.Z A.:. leng) dataF
      processed = case dev of 
                    CPU -> ALN.run1 (((curry pWignerVille) window) . hilbert . supAVG sAVGflag) pData
                    GPU -> ALP.run1 (((curry pWignerVille) window) . hilbert . supAVG sAVGflag) pData
      pList = S.toList $ processed
  file <- openFile newFName WriteMode
  writeString file pList leng 1
  hClose file
  putStrLn "Done !"

makeWVall :: Opts -> IO ()
makeWVall (OptsWV path dev nosAVGflag) = do 
  files <- listDirectory path
  setCurrentDirectory path 
  let fFiles = filter (isSuffixOf ".txt") files
      sAVGflag = A.constant $ not nosAVGflag
  mapM_ (makeWV dev sAVGflag) fFiles
       

makeWV :: CalcDev -> A.Exp Bool -> FilePath -> IO ()
makeWV dev sAVGflag file = do 
  putStrLn $ "processing " ++ file ++ "..."
  text <- TI.readFile file
  let parseRes = parseOnly parseFile text
  case parseRes of 
    (Left errStr) -> error $ errStr
    (Right dataF) -> mapM_ (startWV dev file sAVGflag) dataF
  
startWV :: CalcDev -> FilePath -> A.Exp Bool -> (T.Text,[Double]) -> IO ()
startWV dev oldName sAVGflag (file,dataF) = do 
  let newFName = oldName ++ "-" ++ T.unpack file ++ ".txt"
  putStr $ "   Creating file " ++ newFName ++ " ..."
  let leng = length dataF
      pData = A.fromList (A.Z A.:. leng) dataF
      processed = case dev of 
                    CPU -> ALN.run1 (wignerVille . hilbert . supAVG sAVGflag) pData
                    GPU -> ALP.run1 (wignerVille . hilbert . supAVG sAVGflag) pData
      pList = S.toList $  processed
  file <- openFile newFName WriteMode
  writeString file pList leng 1
  hClose file
  putStrLn "Done !"

supAVG :: A.Exp Bool -> A.Acc (Array DIM1 Double) -> A.Acc (Array DIM1 Double)
supAVG flag arr = 
  A.acond flag (A.map (\x -> x - avg) arr) arr
  where leng = A.length arr
        avg = (A.the $ A.sum arr)/(A.fromIntegral leng)