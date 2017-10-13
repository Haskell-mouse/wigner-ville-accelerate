{-# LANGUAGE OverloadedStrings #-}

module ParseFile where 

import Data.Attoparsec.Text as A
import qualified Data.Text as T
import Data.List
import qualified Control.Applicative as CA 

parseHead :: Parser [T.Text]
parseHead = manyTill (T.pack <$> (many1 (notChar ' ') <* char ' ')) (string "\n" CA.<|> string "\r\n")

parseData :: Parser [[Double]]
parseData = transpose <$> sepBy (sepBy A.double (CA.many $ char ' ')) (CA.many (char ' ') <* (string "\n" CA.<|> string "\r\n"))

parseFile :: Parser [(T.Text,[Double])]
parseFile = zip <$> parseHead <*> parseData

