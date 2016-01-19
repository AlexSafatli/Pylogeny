NAME=Pylogeny
PKG=pylogeny
DOCS=docs
FILES=$(shell find ${PKG} -type f -iname "*.py")
VER=$(shell python -c "from ${PKG}.__version__ import VERSION; print VERSION")
DIST=sdist upload
SDOC=sphinx-apidoc
PY=python
.PHONY: check_environment build tests install dist docs clean

build: ${FILES} setup.py
	${PY} setup.py build

install: build
	-pip uninstall ${PKG}
	${PY} setup.py install

dist: build
	${PY} setup.py ${DIST}

docs: ${FILES} check_environment install
	${SDOC} -F -e -H ${NAME} -A "Alex Safatli" -V ${VER} -R ${VER} -o ${DOCS} ${PKG}
	@cd ${DOCS} && cat conf >> conf.py && make html && make latex && make latexpdf
	cp -r ${DOCS}/_build/html/* $(GH_DOC_REPO)
	cp ${DOCS}/_build/latex/*.pdf docs/api.pdf
	@cd $(GH_DOC_REPO) && git add * && git commit -m "Version ${VER} documentation." && git push origin gh-pages

tests:
	${PY} tests/allTests.py

clean:
	-find ${DOCS} -type f -not -name api.pdf -not -name conf -delete
	-find . -type f -name *.pyc -delete
	-rm -r build dist ${PKG}.egg-info

check_environment:
ifndef GH_DOC_REPO
	$(error GH_DOC_REPO is undefined)
endif
