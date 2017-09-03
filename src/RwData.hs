{-# OPTIONS_HADDOCK ignore-exports #-}
{-# LANGUAGE OverloadedStrings #-}

-- | Module with parsers for command line arguments and cmd help.

module RwData (opts) where

import PseudoWigner
import Data.Text as T
import Options.Applicative as O
import Data.Semigroup ((<>))
default (T.Text)

data WVMode = WV | PWV | SPWV 
data Opts = OptsWV FilePath CalcDev FilePath | OptsPWV FilePath CalcDev Int WindowFunc FilePath | OptsSPWV FilePath CalcDev Int WindowFunc Int WindowFunc FilePath
data CalcDev = CPU | GPU
-- | This function uses execParser function to do all work with catching cmd arguments
-- and parse it with parser that has been created by combine of simple parsers. 

opts :: IO Opts
opts = execParser optI

-- | Parse function . Combine parser for cmd args with info about program. 

optI :: ParserInfo Opts
optI = info ((optWV <|> optPWV) <**> helper) ( header "Test router perfomance through Iperf3")

-- | Parse function for cmd args. It just combine parsers for each option. 

optWV :: O.Parser Opts
optWV = OptsWV 
  <$> wData
  <*> gpu_cpu
  <*> oData

optPWV :: O.Parser Opts  
optPWV = OptsPWV
  <$> wData
  <*> gpu_cpu
  <*> twindow
  <*> winfunc
  <*> oData


wData :: O.Parser FilePath
wData = strOption
 (   long "data"
  <> short 'D'
  <> metavar "FILE"
  <> help "File with data for Wigner-Ville distribution."
 )

oData :: O.Parser FilePath
oData = strOption
 (   long "output"
  <> short 'O'
  <> metavar "FILE"
  <> help "Name of file for output data."
 )

-- | Flag to check if we should clean FW rules before testing. 

pseudo :: O.Parser WVMode
pseudo = flag WV PWV
  (   long "pseudo"  
   <> short 'P'
   <> help "Make peudo Wigner-Ville distribution. You shoud add window leng(Odd) and window function."
  )

smooth_pseudo :: O.Parser WVMode
smooth_pseudo = flag WV SPWV
  (   long "smoothed-pseudo"  
   <> short 'S'
   <> help "Make smoothed peudo Wigner-Ville distribution. You shoud add window leng(Odd) and window function."
  )

gpu_cpu :: O.Parser CalcDev
gpu_cpu = flag CPU GPU
  (   long "gpu"  
   <> short 'G'
   <> help "Set GPU as main calculation device."
  ) 

twindow :: O.Parser Int
twindow = option auto
  (   long "twindow"
   <> short 'T'
   <> metavar "Legth"
   <> showDefault
   <> value 1025
   <> help "Length of window in time domain."
  )

winfunc :: O.Parser WindowFunc
winfunc = option auto 
  (   long "winfunc"
   <> short 'F'
   <> metavar "WinFunc"
   <> showDefault
   <> value Rect
   <> help "Select window function."
  )

