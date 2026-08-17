How I would explain session 3 in three minutes

When we estimate one rate, like the share of fish slipping past at night, we
get one best number: detected over total is 0.110. A range matters more, and the
honest way to get one is to ask which other rates explain the data almost as well.
Plot that and a hump peaks at the estimate; drop a line a fixed amount below the
peak and where it cuts the hump are the endpoints, here [0.082, 0.143], with no
symmetry assumed.

There is more than one recipe for that range: a symmetric delta band from a standard
error, the likelihood read off the curve, a bootstrap of the count, and a Bayesian
blend of prior and data. Here the widest the four disagree at either end is 0.0018,
so they land in nearly the same place. They split only when data are thin or the rate hugs a
boundary, and there the delta band misleads first. Both tools ship the bootstrap.
