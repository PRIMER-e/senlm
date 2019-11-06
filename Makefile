## --- Makefile for Github/R package
## ---

## --- Define R package
PACKAGE=senlm

## --- Define github location
GIT=dciwalsh/senlm

## --- Phony targets
.PHONY: doc

## --- Package creation
## > usethis::create_package ("...PATH.../PACKAGE", rstudio=TRUE)

## --- Initialisation
## git init
## git add DESCRIPTION
## git commit -m "First commit."
## git remote add origin git@github.com:USERNAME/PACKAGE.git
## git push -u origin master

## --- Remove file
## Remove file from git webpage
## git filter-branch
## git pull
## git push

## --- IGNORE FILES
## --- .gitignore
## .Rproj.user
## .Rbuildignore
## senlm.Rproj
## Makefile
## code/*
## *~
## --- .Rbuildignore
## ^.*\.Rproj$
## ^\.Rproj\.user$
## ^Makefile$
## ^code$
## ^run$
## ---

## --- DEFAULT
all :
	@ echo "Do what?"


## --- GIT

## --- Add files
add :
	git add -A .

## --- Commit files
commit :
	git commit -m `date +%Y-%m-%d`

## --- List git files
ls :
	git ls-tree --full-tree -r --name-only HEAD

## --- Push
push :
	git push

## --- Pull
pull :
	git pull

## --- Upload code to github
git :
	make add; make commit; make push


## --- PACKAGE

## --- Documentation
doc :
	echo "devtools::document()" | R --quiet --no-save

## --- Vignettes
vig :
	echo "devtools::clean_vignettes(); devtools::build_vignettes();" | R --quiet --no-save

## --- Local installation
##install :
##	echo "setwd('..'); devtools::install('${PACKAGE}', build_vignettes=TRUE)" | R --quiet --no-save
install :
	echo "devtools::install()" | R --quiet --no-save

## --- Build
build :
	make doc; make vig;  make install

## --- Rebuild, upload to git
update :
	make build; make git;

## --- Get : Reinstall from git
get :
	echo "devtools::install_github('${GIT}', build_vignettes=TRUE)" | R --quiet --no-save


## --- CHECK/TEST

## --- Check
check :
	echo "devtools::load_all(); devtools::check()" | R --quiet --no-save

## --- Test : Initialisation
init-test :
	echo "usethis::use_testthat()" | R --quiet --no-save

## --- Run tests
test :
	echo "devtools::test()" | R --quiet --no-save

## --- Build site
## pkgdown::build_site()
