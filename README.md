This repository provides all the code to replicate the PMCE analysis. 

# REQUIREMENTS
The following R libraries are required to perform the analysis.

* [glmnet, see https://cran.r-project.org/web/packages/glmnet/index.html] with the command:
<pre><code>if (!require("glmnet")) install.packages("glmnet")</code></pre>

* [ggplot2, see https://cran.r-project.org/web/packages/ggplot2/index.html] with the command:
<pre><code>if (!require("ggplot2")) install.packages("ggplot2")</code></pre>

* [gtools, see https://cran.r-project.org/web/packages/gtools/index.html] with the command:
<pre><code>if (!require("gtools")) install.packages("gtools")</code></pre>

* [gridExtra, see https://cran.r-project.org/web/packages/gridExtra/index.html] with the command:
<pre><code>if (!require("gridExtra")) install.packages("gridExtra")</code></pre>

* [igraph, see https://cran.r-project.org/web/packages/igraph/index.html] with the command:
<pre><code>if (!require("igraph")) install.packages("igraph")</code></pre>

* [Rmpfr, see https://cran.r-project.org/web/packages/Rmpfr/index.html] with the command:
<pre><code>if (!require("Rmpfr")) install.packages("Rmpfr")</code></pre>

* [survival, see https://cran.r-project.org/web/packages/survival/index.html] with the command:
<pre><code>if (!require("survival")) install.packages("survival")</code></pre>

* [survminer, see https://cran.r-project.org/web/packages/survminer/index.html] with the command:
<pre><code>if (!require("survminer")) install.packages("survminer")</code></pre>

# RUNNING
The R scripts can be executed either by R GUI or from terminal, with the following commands: 

	Rscript 1_data_processing.R

	Rscript 2_extract_relations.R

	Rscript 3_survival_lasso_analysis.R

Please feel free to contact us if you have problems running our scripts at daniele.ramazzotti1@gmail.com 
