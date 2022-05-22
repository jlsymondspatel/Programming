------------------------------
--  Programming in Haskell  --
--  Chapter 5 Solutions     --
--  J.L. Symonds Patel      --
------------------------------
import Data.Char

----Exercise 1
ex1 = sum [x^2 | x <- [1..100]]

---- Exercise 2
grid :: Int -> Int -> [(Int,Int)]
grid a b = [(x,y) | x <- [0..a], y <- [0..b]]

---- Exercise 3
square :: Int -> [(Int,Int)]
square a = [(x,y) | (x,y) <- (grid a a), x /= y]

---- Exercise 4
replicateNew :: Int -> a -> [a]
replicateNew n x = [x | _ <- [1..n]]

---- Exercise 5
pyths :: Int -> [(Int,Int,Int)]
pyths n = [(x,y,z) | x <- [1..n], y <- [1..n], z <- [1..n], x^2 + y^2 == z^2]

---- Exercise 6
factors :: Int -> [Int]
factors n = [x | x <- [1..n], mod n x == 0]

perfects :: Int -> [Int]
perfects n = [x | x <- [1..n],
              sum (take ((length (factors x))-1) (factors x)) == x]

---- Exercise 7
-- original comprehension:
ex7old = [(x,y) | x <- [1,2], y <- [3,4]]
-- new comprehensions:
ex7new = concat [[(x,y) | y <- [3,4]] | x <- [1,2]]  

---- Exercise 8
find :: Eq a => a -> [(a,b)] -> [b]
find k t = [v | (k',v) <- t, k == k']

positionsNew :: Eq a => a -> [a] -> [Int]
positionsNew x xs = find x (zip xs [0..(length xs)])

---- Exercise 9
scalarproduct :: [Int] -> [Int] -> Int
scalarproduct xs ys = sum [x*y | (x,y) <- zip xs ys]

---- Exercise 10
-- encoding and decoding
let2int :: Char -> Int
let2int c | isLower c = ord c - ord 'a'
          | otherwise = ord c - ord 'A'

int2letlow :: Int -> Char
int2letlow n = chr (ord 'a' + n)

int2letup :: Int -> Char
int2letup n = chr (ord 'A' + n)

shift :: Int -> Char -> Char
shift n c | isLower c = int2letlow (mod (let2int c + n) 26)
          | isUpper c = int2letup (mod (let2int c + n) 26)
          | otherwise = c

encode n xs = [shift n x | x <- xs]

-- frequency tables
table :: [Float]
table = [8.1, 1.5, 2.8, 4.2, 12.7, 2.2, 2.0, 6.1, 7.0,
         0.2, 0.8, 4.0, 2.4, 6.7, 7.5, 1.9, 0.1, 6.0,
         6.3, 9.0, 2.8, 1.0, 2.4, 0.2, 2.0, 0.1]

percent :: Int -> Int -> Float
percent n m = (fromIntegral n / fromIntegral m) * 100

count :: Char -> String -> Int
count x xs = length [x' | x' <- xs, x == x']

letters :: String -> Int
letters xs = length [x | x <- xs,
                     (x >= 'a') && (x <= 'z') || (x >= 'A') && (x <= 'Z')]

freqs :: String -> [Float]
freqs zs = [percent ((count x zs)+(count y zs)) n
           | (x,y) <- zip ['A'..'Z'] ['a'..'z']]
           where n = letters zs

chisqr :: [Float] -> [Float] -> Float
chisqr os es = sum [((o-e)^2)/e | (o,e) <- zip os es]

rotate :: Int -> [a] -> [a]
rotate n xs = drop n xs ++ take n xs

positions :: Eq a => a -> [a] -> [Int]
positions x xs = [i | (x',i) <- zip xs [0..], x == x']

crack :: String -> String
crack xs = encode (-factor) xs
           where
              factor = head (positions (minimum chitab) chitab)
              chitab = [chisqr (rotate n table') table | n <- [0..25]]
              table' = freqs xs
