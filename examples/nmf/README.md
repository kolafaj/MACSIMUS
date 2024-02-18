# Normal mode vibrations of biphenyl

## Copy biphenyl (or whatever you wish)

From this directory:

`$ cp ../../blend/che/biphenyl.che .`

Alternatively, you can work in `blend/che/`

`$ cd ../../blend/che/`

## Run script

Three options:

* `$ nmf.sh biphenyl.che`
* `doubleclick biphenyl.che` from Midnight Commander and select 1
* `start biphenyl.che` and select 1

Then, check that the molecule is optimized:

* Click [CG] or type [,] in panel `--- optimize ---`
* Optionally, push [ran] or type ':'
* Click [finish] or type [.]

A window of `show` will appear:

* Use buttons [<] and [>] or keys PgUp and PgDn to switch modes
* Modes 0..5 (or 0..4 for a linear molecule) correspond to translation/rotation (zero eigenvalues)
