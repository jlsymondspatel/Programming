------------------------------
--  Programming in Haskell  --
--  Chapter 4 Solutions     --
--  J.L. Symonds Patel      --
------------------------------

-- Exercise One
halve :: [a] -> ([a],[a])
halve xs | mod (length xs) (2) == 0 = (take (div (length xs) 2) xs,drop (div (length xs) 2) xs)
         | otherwise = ([],[])

-- Exercise Two
thirdht :: [a] -> a
thirdht xs = head (tail (tail xs))

thirdind :: [a] -> a
thirdind xs = xs !! 2

thirdpat :: [a] -> a
thirdpat (_:_:x:_) = x

-- Exercise Three
safetailif :: [a] -> [a]
safetailif xs = if null xs
                then []
                else tail xs

safetailguard :: [a] -> [a]
safetailguard xs | not (null xs) = tail xs
                 | otherwise = []

safetailpat :: [a] -> [a]
safetailpat [] = []
safetailpat xs = tail xs

-- Exercise Four
disj1 :: Bool -> Bool -> Bool
disj1 True True = True
disj1 True False = True
disj1 False True = True
disj1 False False = False

disj2 :: Bool -> Bool -> Bool
disj2 False a = a
disj2 True _ = True

disj3 :: Bool -> Bool -> Bool
disj3 False False = False
disj3 _ _ = True

disj4 :: Bool -> Bool -> Bool
disj4 a b | a == b = a
          | otherwise = False

-- Exercise Five
(&&*) :: Bool -> Bool -> Bool
a &&* b = if a == True
           then if b == True
                then True
                else False
           else False

-- Exercise Six
(&&**) :: Bool -> Bool -> Bool
a &&** b = if a == b
          then a
          else False

-- Exercise Seven
mult :: Int -> (Int -> (Int -> Int))
mult = \x -> (\y -> (\z -> x*y*z))

-- Exercise Eight
luhnDouble :: Int -> Int
luhnDouble x | (2*x) > 9 = (2*x)-9
             | otherwise = 2*x

luhn :: Int -> Int -> Int -> Int -> Bool
luhn w x y z = if mod ((luhnDouble w) + x + (luhnDouble y) + z) 10 == 0
               then True
               else False
