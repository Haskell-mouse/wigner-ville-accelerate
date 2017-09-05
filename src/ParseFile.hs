{-# LANGUAGE OverloadedStrings #-}

module ParseFile where 

import Data.Attoparsec.Text as A
import qualified Data.Text as T
import Data.List

parseHead :: Parser [T.Text]
parseHead = manyTill (T.pack <$> (many1 (notChar ' ') <* char ' ')) (char '\n')

parseData :: Parser [[Double]]
parseData = transpose <$> sepBy (sepBy A.double (char ' ')) (string " \n")

parseFile :: Parser [(T.Text,[Double])]
parseFile = zip <$> parseHead <*> parseData

