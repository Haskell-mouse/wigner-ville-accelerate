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
import Data.Monoid
import Debug.Trace 
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
       (OptsWV _ _) -> makeWVall args
       (OptsPWV _ _ _ _) -> makePWVall args
  
writeString :: Handle -> [Double] -> Int -> Int -> IO ()
writeString _ [] _ _  = return ()
writeString file (x:xs) w n = do 
  TI.hPutStr file $ DT.toPrecision 5 x <> " "
  if (n `mod` w) == 0 
  then do 
    TI.hPutStr file "\n"
    writeString file xs w (n + 1)
  else writeString file xs w (n + 1)

makePWVall :: Opts -> IO ()
makePWVall (OptsPWV path dev wlen wfun) = do
  if (even wlen)
  then error "Window length maust be odd !!"
  else do  
    files <- listDirectory path
    setCurrentDirectory path 
    let fFiles = filter (isSuffixOf ".txt") files
    let window = makeWindow wfun (A.unit $ A.constant wlen) :: A.Acc (A.Array A.DIM1 Double)
    mapM_ (makePWV dev window) fFiles

makePWV :: CalcDev -> A.Acc (A.Array A.DIM1 Double) -> FilePath -> IO ()
makePWV dev window file = do 
  putStrLn $ "processing " ++ file ++ "..."
  text <- TI.readFile file
  let (Right dataF) = parseOnly parseFile text
  mapM_ (startPWV dev window file) dataF

startPWV :: CalcDev -> A.Acc (A.Array A.DIM1 Double) -> FilePath -> (T.Text,[Double]) -> IO ()
startPWV dev window oldName (file,dataF) = do
  let newFName = oldName ++ "-" ++ T.unpack file ++ ".txt" ++ " ..."
  putStr $ "   Creating file " ++ newFName 
  let leng = length dataF
      pData = A.fromList (A.Z A.:. leng) dataF
      processed = case dev of 
                    CPU -> ALN.run1 (((curry pWignerVille) window) . hilbert) pData
                    GPU -> ALP.run1 (((curry pWignerVille) window) . hilbert) pData
      pList = S.toList processed
  file <- openFile newFName WriteMode
  writeString file pList leng 1
  hClose file
  putStrLn "Done !"

makeWVall :: Opts -> IO ()
makeWVall (OptsWV path dev) = do 
  files <- listDirectory path
  setCurrentDirectory path 
  let fFiles = filter (isSuffixOf ".txt") files
  mapM_ (makeWV dev) fFiles
       

makeWV :: CalcDev -> FilePath -> IO ()
makeWV dev file = do 
  putStrLn $ "processing " ++ file ++ "..."
  text <- TI.readFile file
  let (Right dataF) = parseOnly parseFile text
  mapM_ (startWV dev file) dataF
  
startWV :: CalcDev -> FilePath -> (T.Text,[Double]) -> IO ()
startWV dev oldName (file,dataF) = do 
  let newFName = oldName ++ "-" ++ T.unpack file ++ ".txt"
  putStr $ "   Creating file " ++ newFName ++ " ..."
  let leng = length dataF
      pData = A.fromList (A.Z A.:. leng) dataF
      processed = case dev of 
                    CPU -> ALN.run1 (wignerVille . hilbert) pData
                    GPU -> ALP.run1 (wignerVille . hilbert) pData
      pList = S.toList processed
  file <- openFile newFName WriteMode
  writeString file pList leng 1
  hClose file
  putStrLn "Done !"