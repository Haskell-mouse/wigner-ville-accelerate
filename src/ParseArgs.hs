{-# OPTIONS_HADDOCK ignore-exports #-}
{-# LANGUAGE OverloadedStrings #-}

-- | Module with parsers for command line arguments and cmd help.

module ParseArgs(WVMode(..), opts, Opts(..), CalcDev(..)) where

import PseudoWigner
import Data.Text as T
import Options.Applicative as O
import Data.Semigroup ((<>))
default (T.Text)

data WVMode = WV | PWV | SPWV 
type NoSubAVG = Bool
data Opts = OptsWV FilePath CalcDev NoSubAVG | OptsPWV FilePath CalcDev Int WindowFunc NoSubAVG | OptsSPWV FilePath CalcDev Int WindowFunc Int WindowFunc NoSubAVG
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
  <$> dataPath
  <*> gpu_cpu
  <*> no_subtract_avg

optPWV :: O.Parser Opts  
optPWV = pseudo 
   *> (OptsPWV
  <$> dataPath
  <*> gpu_cpu
  <*> twindow
  <*> winfunc
  <*> no_subtract_avg )

dataPath :: O.Parser FilePath
dataPath = strOption
 (   long "data"
  <> short 'D'
  <> metavar "PATH"
  <> help "Path with files for Wigner-Ville distribution."
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
   <> help "Length of frequency smoothing window in time domain."
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

no_subtract_avg :: O.Parser NoSubAVG
no_subtract_avg = switch
  (   long "no-subtract-average"
   <> short 'N'
   <> help "If enabled, then dont subtract average value from all elemets of array before Hilbert transform."
  )