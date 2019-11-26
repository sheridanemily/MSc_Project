# README

This subdirectory contains the example scripts to accompany the manuscript.
Read the section [Version your code][version] for the full explanation.
The commands described below should be executed from within the subdirectory `code`.

[version]: https://www.authorea.com/users/5990/articles/17489/_show_article#article-paragraph-version__minus__your__minus__code__dot__tex-landing-welcome

## process.sh

`process.sh` downloads the ENCODE CTCF ChIP-seq data from multiple types of kidney samples and calls peaks.
The downloaded BAM files and [MACS2][] output are saved in tissue-specific subdirectories in `../data/`.
Run it as follows:

```bash
bash process.sh
```

You will need to have Python 2.7 installed because it is a [prerequisite][] for the peak caller [MACS2][].
To install MACS2, run:

```bash
pip install macs2
```

[prerequisite]: https://github.com/taoliu/MACS/blob/master/INSTALL.rst#prerequisites
[MACS2]: https://github.com/taoliu/MACS

## clean.py

`clean.py` filters peaks with a fold change cutoff and merges peaks from the different kidney samples.
It creates the file `../data/sites-union.bed`.
Run it as follows:

```bash
python clean.py
```

It can be run with either Python2 or Python3.
You will first need to install [bedtools][] and [pybedtools][].
See their documentation for instructions, but as an example, this can be accomplished on Ubuntu with the following:

```bash
apt-get install bedtools
pip install pybedtools
```

[bedtools]: http://bedtools.readthedocs.org/en/latest/content/installation.html
[pybedtools]: http://pythonhosted.org/pybedtools/main.html

## analyze.R

`analyze.R` creates diagnostic plots on the length of the peaks and their distribution across the genome.
It creates the file `../data/sites-union.pdf`.
Run it as follows:

```bash
Rscript analyze.R
```

sudo apt-get update
sudo apt-get install git -y

#### Installing git

git config --global user.name "sheridanemily"
git config --global user.email "sheridanemilyanne@gmail.com"

#### Providing username and email of user

cd /Downloads
tar xvf codes.tar
mkdir -p ../thesis
cp code/* ../thesis
cd ../thesis
ls

#### Download tar files and store in appropriate file system, in this case I store them in the /thesis directory.

git init
git config user.name "sheridanemily"
git config user.email "sheridanemilyanne@gmail.com"
git status

#### Configure git locally and check the staus

git add process.sh
git status
git add clean.py analyze.R
git -commit -m "Add initial version of thesis code."
git log

#### Stage and commit files to the local repo

(base) a@e:~/thesis$ git log
commit a55ee6ac0095877739fc23610c37b0817e9f9f7c (HEAD -> master)
Author: sheridanemily <sheridanemilyanne@gmail.com>
Date:   Tue Nov 26 13:38:48 2019 +0000

    Add initial version of thesis code.
    
#### The 'git log' function provides basic information about the log. It shows who is the author of the log, when it was created, and how many changes have been made to the log

nano -w clean.py

####### Change fc_cutoff in clean.py
#+fc_cutoff = 20
#-fc_cutoff = 10

#### This function allows the author to make changes to the clean.py file

(base) a@e:~/thesis$ git diff
diff --git a/clean.py b/clean.py
index 9f76681..20384b2 100644
--- a/clean.py
+++ b/clean.py
@@ -28,7 +28,8 @@ def filter_fold_change(feature, fc = 1):
         return False
 
 #Filter based on fold-change over control sample
-fc_cutoff = 10
+fc_cutoff = 20
 epithelial = epithelial.filter(filter_fold_change, fc = fc_cutoff).saveas()
 proximal_tube = proximal_tube.filter(filter_fold_change, fc = fc_cutoff).saveas()
 kidney = kidney.filter(filter_fold_change, fc = fc_cutoff).saveas()

#### The purpose of the 'git diff' function is to view the changes you made relative to the index. In this case, the original fc_cutoff was 10 and it was changed to 20


https://github.com/sheridanemily/thesis.git


