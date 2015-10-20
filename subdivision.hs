{-|
Module      : Subdivision
Description : Univariate subdivisions
Exercises from Analysis and Design of Univariate Subdivision Schemes (by M. Sabin).
-}
module Subdivision where

import Data.List (tails, transpose)
import Data.List.Split (chunksOf)

-- * Part II. Dramatis Personae
-- ** Section 11: An introduction to some regularly-appearing characters

-- | A subdivision scheme is defined by its mask and its arity.
-- The mask can be given in integers, as in @Scheme 2 [1,3,3,1]@.
-- A suitable divisor for normalizing is computed by 'divisor'.
data Scheme = Scheme { arity :: Int, mask  :: [Double] }
              deriving (Show)

-- | Cubic B-spline: @[1,4,6,4,1]/8@
cubicBSpline            = Scheme 2 [1,4,6,4,1]
-- | Quadratic B-spline: @[1,3,3,1]/4@
quadraticBSline         = Scheme 2 [1,3,3,1]
-- | Ternary Quadratic B-spline: @[1,3,6,7,6,3,1]/9@
ternaryQuadraticBSpline = Scheme 3 [1,3,6,7,6,3,1]
-- | Ternary neither: @[1,3,5,5,3,1]/6@
ternaryNeither          = Scheme 3 [1,3,5,5,3,1]
-- | Four-point: @[-1,0,9,16,9,0,-1]/16@
fourPoint               = Scheme 2 [-1,0,9,16,9,0,-1]

-- | A list of the above schemes, for convenience.
schemes = [ cubicBSpline
          , quadraticBSline
          , ternaryQuadraticBSpline
          , ternaryNeither
          , fourPoint
          ]

-- | @'divisor' s@ Returns the number with which the mask is divided.
--
-- >>> divisor cubicBSpline
-- 8.0
--
-- prop> divisor s = sum (mask s) / fromIntegral (arity s)
divisor :: Scheme -> Double
divisor (Scheme a xs) = sum xs / fromIntegral a

-- | This is the mask to be used in computations.
dividedMask :: Scheme -> [Double]
dividedMask s@(Scheme a xs) = map (/ d) xs
    where d = divisor s

-- | The result is a list of @a@ lists, where @a@ is the arity of the scheme.
stencils :: Scheme -> [[Double]]
stencils (Scheme a xs) = transpose (chunksOf a xs)

-- * Part III. Analyses
-- ** Section 12: Support

-- | Uses the scheme on itself (interpreted as 1-dimensional control points).
-- This results in a more dense scheme with the same support.
-- Iterated use of this function gives a good approximation of the basis function.
square :: Scheme -> Scheme
square (Scheme a xs) = Scheme (a * a) (foldr f start xs)
    where n      = length xs
          start  = replicate (n - a) 0
          f x ys = let (xs1,xs2) = splitAt a $ map (* x) xs
                       (ys1,ys2) = splitAt (n - a) ys
                   in xs1 ++ (zipWith (+) xs2 ys1) ++ ys2

-- | Returns the support of the scheme. Note that it may not be integral:
--
-- >>> support ternaryNeither
-- 2.5
support :: Scheme -> Double
support (Scheme a xs) = fromIntegral (length xs - 1) / fromIntegral (a - 1)

-- | @'practicalSupports' n s@ approximates the support of scheme @s@,
-- where the basis function has considerable impact.
-- The result is a list of three values corresponding to the tolerances 1%, 2%, and 5%:
--
-- >practicalSupports 3 cubicBSpline == [3.21484375,3.01171875,2.66015625]
--
-- It uses @n@ squaring operations to generate a good approximation of the basis function.
practicalSupports :: Int -> Scheme -> [Double]
practicalSupports n = largeRatios . (!! n) . iterate square
    where largeRatios s     = map (largeRatio s) [0.01, 0.02, 0.05]
          largeRatio  s tol = fromIntegral (largeElements s) / fromIntegral (arity s)
              where largeElements = length . filter (> tol) . dividedMask

-- ** Section 13: Enclosure

-- | @'norm' n s@ approximates the l-infinity norm of the scheme @s@.
-- It uses @n@ squaring operations to generate a good approximation of the basis function.
--
-- >>> norm 3 fourPoint
-- 1.2511177496053278
-- >>> norm 4 (Scheme 2 [1,2,1,-1,-2,-1,1,2,1])
-- 1.0001220703125
norm :: Int -> Scheme -> Double
norm n = maxRow . rows . (!! n) . iterate square
    where rows s = transpose . chunksOf (arity s) . dividedMask $ s
          maxRow = maximum . map (sum . map abs)
