{-# LANGUAGE BangPatterns #-}
module RwData where

import qualified Data.Text.IO as TI
import qualified Data.Text as T
import qualified Data.Vector.Storable  as U
import qualified Data.Array.Accelerate as A
import qualified Data.Array.Accelerate.IO as AI
import Data.Either
import Data.Char 
import qualified Data.Text.Read as R
import System.Environment

type Period = Double

parse :: [String] ->IO (Maybe ([A.Array A.DIM1 Double]))
parse args 
    | args == [] = do
        putStrLn "usage: wigner [optional : '-p winfile(or winlength)' or '-sp winfile1(or timeWinLength) winfile2(or timeWinLength)'] datafile "
        return Nothing
    | (head args == "-h") || (head args == "--help") = do
	    putStrLn "usage: wigner [optional : '-p winfile(or winlength)' or '-sp winfile1(or timeWinLength) winfile2(or timeWinLength)'] datafile "
	    return Nothing
    | (length args == 1) && (head args /= "-p") && (head args /= "-sp") = do
    	dataV <- prepareDataWV (args)
        return $! Just ([dataV])
    | (head args == "-p") && (length args == 3) = do
    	allData <- prepareDataPWV (init args)
        return $! Just (allData)
    | (head args == "-sp") && (length args == 4) = do
        allData <- prepareDataSPWV args
        return $! Just (allData)
    | otherwise = do
	    putStrLn "usage: wigner [optional : '-p winfile(or winlength)' or '-sp winfile1(or timeWinLength) winfile2(or timeWinLength)']  datafile "
	    return Nothing

prepareDataWV :: [String] -> IO (A.Array A.DIM1 Double)
prepareDataWV args = do
     let filename = head args
     dataV <- readVectFile filename
     return $! dataV


prepareDataPWV :: [String] -> IO [A.Array A.DIM1 Double]
prepareDataPWV args = do
    let timeWindowFile = head args
    let dataFile = head (tail args)
    dataV <- readVectFile dataFile
    tWindow <- readVectFile timeWindowFile
    return $! [tWindow,dataV]

prepareDataSPWV :: [String] -> IO [A.Array A.DIM1 Double]
prepareDataSPWV args = do
    allData <- mapM readVectFile $ tail args
    return $! allData

readVectFile :: String -> IO (A.Array A.DIM1 Double)
readVectFile file = do
	s <- TI.readFile file
	let vect = parseText (T.snoc s ' ')
	let vlength = U.length vect 
	return $ AI.fromVectors (A.Z A.:. vlength) vect 

parseText :: T.Text -> U.Vector Double
parseText = U.unfoldr step
  where
     step !s = case R.double s of
        Right (!k, !t) -> Just (k,T.tail t)
        Left _ -> Nothing 

newVect :: Int -> A.Array A.DIM1 Double
newVect num 
    | even num = let list = take num $ repeat 1
                     lengthL = length list
                 in  A.fromList (A.Z A.:. lengthL) list
	| odd num  = newVect (num+1)