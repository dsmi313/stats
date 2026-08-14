How I would explain session 7 in three minutes

When a number we need, like guidance efficiency, is uncertain, we do not want a
single plugged-in value; we want to account for every value it might take,
weighted by how likely each is. Written as math that is an integral, and it can
look intimidating.

There is a shortcut that gives the identical answer. Draw a few thousand values
of the uncertain number from its range, run the escapement calculation once for
each, and look at the spread of results. The average lands on the integral and
the middle band is the honest interval. That is all the bootstrap loop is doing
when it reads one guidance-efficiency draw and one stock draw per pass: it is
integrating out our uncertainty by drawing, not by calculus.
