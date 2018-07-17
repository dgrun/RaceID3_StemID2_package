# RaceID algorithm

RaceID is a clustering algorithm for the identification of cell types
from single-cell RNA-sequencing data. It was specifically designed for
the detection of rare cells which correspond to outliers in
conventional clustering methods. The package contains RaceID3, the
most recently published version of this algorithm, and StemID2, an
algorithm for the identification of lineage trees based on RaceID3
analysis. RaceID3 utilizes single cell expression data, and was
designed to work well with quantitative single-cell RNA-seq data
incorporating unique molecular identifiers. It requires a gene-by-cell
expression matrix as input and produces a clustering partition
representing cell types. StemID2 assembles these cell types into a
lineage tree.

## Installing

After downloading and unzipping
```
unzip RaceID3_StemID2_package-master.zip 
```

it can be installed from the command line by
```
R CMD INSTALL RaceID3_StemID2_package-master
```

or directly in R from source by
```
install.packages("RaceID3_StemID2_package-master",repos = NULL, type="source")
```
(if R is started from the directory where `RaceID3_StemID2_package-master.zip` has been downloaded to; otherwise specify the full path)


Alternatively, install in R directly from github using devtools:
```
install.packages("devtools")
library(devtools)
install_github("dgrun/RaceID3_StemID2_package")
```

## Running a RaceID analysis

Load package:
```
library(RaceID)
```

See vignette for details and examples:
```
vignette("RaceID")
```

## Reference:

Herman JS, Sagar, Grün D. (2018) FateID infers cell fate bias in
multipotent progenitors from single-cell RNA-seq data. Nat Methods. 2018 May;15(5):379-386. doi: 10.1038/nmeth.4662.
