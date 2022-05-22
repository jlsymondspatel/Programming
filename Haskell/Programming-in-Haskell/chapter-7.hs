------------------------------
--  Programming in Haskell  --
--  Chapter 7 Solutions     --
--  J.L. Symonds Patel      --
------------------------------
import Data.Char

---- Exercise 1:
-- original list comprehension:
-- [f x | x <- xs, p x]
applyFunctionToParticularElementsOfList :: (a -> a) -> (a -> Bool) -> [a] -> [a]
applyFunctionToParticularElementsOfList f p xs = map f (filter p xs)

---- Exercise 2:
-- a
all' :: (a -> Bool) -> [a] -> Bool
all' p = and . map p

-- b
any' :: (a -> Bool) -> [a] -> Bool
any' p = or . map p

-- c
takeWhile' :: (a -> Bool) -> [a] -> [a]
takeWhile' _ [] = []
takeWhile' p (x:xs) | p x = x : takeWhile' p xs
                    | otherwise = []

-- d
dropWhile' :: (a -> Bool) -> [a] -> [a]
dropWhile' _ [] = []
dropWhile' p (x:xs) | p x = dropWhile' p xs
                    | otherwise = []

---- Exercise 3
-- map f
map' :: (a -> b) -> [a] -> [b]
map' f = foldr (\x xs -> f x : xs) []

-- filter p
filter' :: (b -> Bool) -> [b] -> [b]
filter' p = foldr (\x xs -> if p x then x : xs else xs) []

---- Exercise 4
dec2int :: [Int] -> Int
dec2int = foldl (\x y -> (10*x)+y) 0

---- Exercise 5
curry' :: ((a,b) -> c) -> (a -> b -> c)
curry' f = (\x y -> f (x,y))

uncurry' :: (a -> b -> c) -> ((a,b) -> c)
uncurry' f = (\(x,y) -> f x y)

---- Exercise 6
--unfold :: (a -> Bool) -> (a -> a) -> (a -> a) -> a -> [a]
-- ^^^ this doesn't work, type mismatch is inevitable.
unfold p h t x | p x = []
               | otherwise = h x : unfold p h t (t x)
-- int2bin
int2bin' :: Int -> [Int]
int2bin' = unfold (== 0) (\x -> mod x 2) (\x -> div x 2)

-- chop8
chop8' :: [Int] -> [[Int]]
chop8' = unfold (null) (take 8) (drop 8)

-- map f
map'' :: (a -> b) -> ([a] -> [b])
map'' f = unfold (null) (f . head) (tail) 

-- iterate
iterate' :: (a -> a) -> a -> [a]
iterate' = unfold (const False) (id)

---- Exercise 7
-- ORIGINAL TRANSMITTER CODE PLUS MODIFICATION
-- Binary string transmitter example from chapter 7 of Programming
-- in Haskell, Graham Hutton, Cambridge University Press, 2016.

-- Base conversion

type Bit = Int

bin2int :: [Bit] -> Int
bin2int = foldr (\x y -> x + 2*y) 0

int2bin :: Int -> [Bit]
int2bin 0 = []
int2bin n = n `mod` 2 : int2bin (n `div` 2)

make8 :: [Bit] -> [Bit]
make8 bits = take 8 (bits ++ repeat 0)

addParity :: [Bit] -> [Bit]
addParity bits | odd (length (filter (== 1) bits)) = [1] ++ bits
               | otherwise = [0] ++ bits

parityUnchanged :: [Bit] -> Bool
parityUnchanged bits | odd (length (filter (== 1) (tail bits))) && (1 == head bits) = True
                     | otherwise = False

parityFilter :: [[Bit]] -> [[Bit]]
parityFilter bitss = if and [parityUnchanged bits | bits <- bitss]
                     then [tail bits | bits <- bitss]
                     else error "parity mismatch"

takeParity :: [Bit] -> [Bit]
takeParity bits = tail bits

-- Transmission

encode :: String -> [Bit]
encode = concat . map (addParity . make8 . int2bin . ord)

chop9 :: [Bit] -> [[Bit]]
chop9 []   = []
chop9 bits = take 9 bits : chop9 (drop 9 bits)

decode :: [Bit] -> String
decode = map (chr . bin2int) . parityFilter . chop9

transmit :: String -> String
transmit = decode . channel . encode

channel :: [Bit] -> [Bit]
channel = id

---- Exercise 8
faultyChannel :: [Bit] -> [Bit]
faultyChannel = tail

faultyTransmit :: String -> String
faultyTransmit = decode . faultyChannel . encode

---- Exercise 9
altMap :: (a -> b) -> (a -> b) -> [a] -> [b]
altMap _ _ [] = []
altMap f g (x:xs) = f x : altMap g f xs

---- Exercise 10
luhnDouble :: Int -> Int
luhnDouble x | (2*x) > 9 = (2*x)-9
             | otherwise = 2*x

applyLuhnDouble :: [Int] -> [Int]
applyLuhnDouble = altMap (id) (luhnDouble) . reverse

luhn :: [Int] -> Bool
luhn xs = if (sum (applyLuhnDouble xs) `mod` 10) == 0
          then True
          else False
