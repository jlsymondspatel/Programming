------------------------------
--  Programming in Haskell  --
--  Chapter 8 Solutions     --
--  J.L. Symonds Patel      --
------------------------------

---- Exercise 1
data Nat = Zero | Succ Nat
         deriving Show

add :: Nat -> Nat -> Nat
add Zero n = n
add (Succ m) n = Succ (add m n)

mult :: Nat -> Nat -> Nat
mult Zero m = Zero
mult (Succ m) n = add (mult m n) n

---- Exercise 2
data Tree a = Leaf a | Node (Tree a) a (Tree a)
            deriving Show

occurs :: Ord a => a -> Tree a -> Bool
occurs x (Leaf y) = x == y
occurs x (Node l y r) = case (compare x y) of
                          LT -> occurs x l
                          EQ -> True
                          GT -> occurs x r

-- This new version requires only one comparison evaluation, as opposed to more.

---- Exercise 3
data Tree2 a = Leaf2 a | Node2 (Tree2 a) (Tree2 a)
             deriving Show

leaves :: Tree2 a -> Int
leaves (Leaf2 _) = 1
leaves (Node2 l r) = (leaves l) + (leaves r)

balanced :: Tree2 a -> Bool
balanced (Leaf2 _) = True
balanced (Node2 l r) = (abs ((leaves l) - (leaves r)) <= 1)
                       && (balanced l) && (balanced r)

---- Exercise 4

halve :: [a] -> ([a],[a])
halve xs = splitAt (div (length xs) 2) xs

balance :: [a] -> Tree2 a
balance [x] = (Leaf2 x)
balance xs = Node2 (balance l) (balance r)
             where l = fst (halve xs)
                   r = snd (halve xs)

---- Exercise 5
data Expr = Val Int | Add Expr Expr

folde :: (Int -> a) -> (a -> a -> a) -> Expr -> a
folde f g (Val x) = f x
folde f g (Add x y) = g (folde f g x) (folde f g y)

---- Exercise 6
eval :: Expr -> Int
eval = folde id (+)

size :: Expr -> Int
size = folde (const 1) (+)

---- Exercise 7
{-
instance Eq a => Eq (Maybe a) where
  Nothing == Nothing = True
  (Just x == Just y) = (x == y)
  _ == _ = False
-}

{-
instance Eq a => Eq [a] where
  (x:xs) == (y:xs) = (x == y) && (xs == ys)
  [] == [] = True
  _ == _ = False
-}

---- Exercise 8
-- Original script plus extensions:
-- Caeser cipher example from chapter 8 of Programming in Haskell,
-- Graham Hutton, Cambridge University Press, 2016.

-- Propositions

data Prop = Const Bool
          | Var Char
          | Not Prop
          | And Prop Prop
          | Imply Prop Prop
          | Or Prop Prop
          | Equiv Prop Prop

-- Substitutions

type Subst = Assoc Char Bool

type Assoc k v = [(k,v)]

find :: Eq k => k -> Assoc k v -> v
find k t = head [v | (k',v) <- t, k == k']

-- Tautology checker

eval2 :: Subst -> Prop -> Bool
eval2 _ (Const b)   = b
eval2 s (Var x)     = find x s
eval2 s (Not p)     = not (eval2 s p)
eval2 s (And p q)   = eval2 s p && eval2 s q
eval2 s (Imply p q) = eval2 s p <= eval2 s q
eval2 s (Or p q)    = eval2 s p || eval2 s q
eval2 s (Equiv p q) = eval2 s p == eval2 s q

vars :: Prop -> [Char]
vars (Const _)   = []
vars (Var x)     = [x]
vars (Not p)     = vars p
vars (And p q)   = vars p ++ vars q
vars (Imply p q) = vars p ++ vars q
vars (Or p q)    = vars p ++ vars q
vars (Equiv p q) = vars p ++ vars q

bools :: Int -> [[Bool]]
bools 0 = [[]]
bools n = map (False:) bss ++ map (True:) bss
          where bss = bools (n-1)

rmdups :: Eq a => [a] -> [a]
rmdups []     = []
rmdups (x:xs) = x : filter (/= x) (rmdups xs)

substs :: Prop -> [Subst]
substs p = map (zip vs) (bools (length vs))
           where vs = rmdups (vars p)

isTaut :: Prop -> Bool
isTaut p = and [eval2 s p | s <- substs p]

---- Exercie 9
-- Original script plus extensions:
-- Abstract machine example from chapter 8 of Programming in Haskell,
-- Graham Hutton, Cambridge University Press, 2016.

data Expr2 = Val2 Int | Add2 Expr2 Expr2 | Mul2 Expr2 Expr2

type Cont = [Op]

data Op = EVAL_ADD Expr2 | EVAL_MUL Expr2 | ADD Int | MUL Int

eval3 :: Expr2 -> Cont -> Int
eval3 (Val2 n)   c = exec c n
eval3 (Add2 x y) c = eval3 x (EVAL_ADD y : c)
eval3 (Mul2 x y) c = eval3 x (EVAL_MUL y : c)

exec :: Cont -> Int -> Int
exec []           n = n
exec (EVAL_ADD y : c) n = eval3 y (ADD n : c)
exec (EVAL_MUL y : c) n = eval3 y (MUL n : c)
exec (ADD n : c) m = exec c (n+m)
exec (MUL n : c) m = exec c (n*m)

value :: Expr2 -> Int
value e = eval3 e []
