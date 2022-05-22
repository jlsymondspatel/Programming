------------------------------
--  Programming in Haskell  --
--  Chapter 6 Solutions     --
--  J.L. Symonds Patel      --
------------------------------

---- Exercise 1:
-- original recursive factorial function
fac0 :: Int -> Int
fac0 0 = 1
fac0 n = n * fac0 (n-1)
-- If a negative number is entered as input, the function continues endlessly since the base case is never met.
-- new factorial function prohibiting negative arguments
fac1 :: Int -> Int
fac1 0 = 1
fac1 n | n > 0 = n * fac0 (n-1)
       | otherwise = 0

---- Exercise 2:
sumdown :: Int  -> Int
sumdown 0 = 0
sumdown n | n > 0 = n + sumdown (n-1)
          | otherwise = 0

---- Exercise 3:
(^*) :: Int -> Int -> Int
_ ^* 0 = 1
m ^* n = m * (m ^ (n-1))
-- how 2 ^* is evaluated
-- 2 ^* 3 = 2 * (2 ^* 2)
--        = 2 * (2 * (2 ^* 1))
--        = 2 * (2 * (2 * (2 ^* 0)))
--        = 2 * (2 * (2 * (1)))

---- Exercise 4:
euclid :: Int -> Int -> Int
euclid x y | x == y = x
           | x < y = euclid x (y-x)
           | x > y = euclid (x-y) y
           | (x < 0) || (y < 0) = 0

---- Exercise 5:
-- length [1,2,3]
{-
= 1 + length [2,3]
= 1 + 1 + length [3]
= 1 + 1 + 1 + length []
= 1 + 1 + 1 + 0
-}

-- drop 3 [1,2,3,4,5]
{-
= drop 2 [2,3,4,5]
= drop 1 [3,4,5]
= drop 0 [4,5]
= [4,5]
-}

-- init [1,2,3]
{-
= 1 : init [2,3]
= 1 : 2 : init [3]
= 1 : 2 : []
= [1,2]
-}

---- Exercise 6:
-- a
and0 :: [Bool] -> Bool
and0 [] = True
and0 (x:xs) = x && (and0 xs)

-- b
concat0 :: [[a]] -> [a]
concat0 [[]] = []
concat0 (x:xs) = x ++ concat0 xs

-- c
replicate0 :: Int -> a -> [a]
replicate0 0 _ = []
replicate0 n m = m : replicate0 (n-1) m

-- d
(!!*) :: [a] -> Int -> a
(x:xs) !!* 0 = x
(x:xs) !!* n = xs !!* (n-1)

-- e
elem0 :: Eq a => a -> [a] -> Bool
elem0 _ [] = False
elem0 n (x:xs) = (x == n) || (elem0 n xs)

---- Exercise 7:
merge :: Ord a => [a] -> [a] -> [a]
merge xs [] = xs
merge [] xs = xs
merge (x:xs) (y:ys) | x < y = x : merge xs (y:ys)
                    | x > y = y : merge (x:xs) ys
                    | otherwise = x : y : merge xs ys

---- Exercise 8:
halve :: [a] -> ([a],[a])
halve xs = splitAt (div (length xs) 2) xs

msort :: Ord a => [a] -> [a]
msort [] = []
msort [x] = [x]
msort xs = merge (msort left) (msort right)
           where (left,right) = halve xs

---- Exercise 9:
-- a
sum0 :: Num a => [a] -> a
sum0 [] = 0
sum0 (x:xs) = x + sum0 xs

-- b
take0 :: Int -> [a] -> [a]
take0 0 _ = []
take0 n (x:xs) = x : take0 (n-1) xs

-- c
last0 :: [a] -> a
last0 [x] = x
last0 (x:xs) = last0 xs
