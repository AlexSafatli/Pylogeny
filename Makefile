NAME=Pylogeny
PKG=pylogeny
DOCS=docs
FILES=$(shell find ${PKG} -type f -iname "*.py")
TESTS=$(shell find tests -type -f -iname "*.py")
VER=$(shell python -c "from ${PKG}.__version__ import VERSION; print VERSION")
DIST=sdist upload
EDOC=epydoc
SDOC=sphinx-apidoc
PY=python
.PHONY: check_environment build install dist docs clean

build: ${FILES} setup.py
	${PY} setup.py build

install: build
	-sudo pip uninstall ${PKG}
	sudo ${PY} setup.py install

dist: build
	${PY} setup.py ${DIST}

docs: ${FILES} check_environment install
	${EDOC} --pdf --name ${NAME} --url "http://AlexSafatli.github.io/Pylogeny" --graph all ${FILES} -o ${DOCS}
	${SDOC} -F -e -H ${NAME} -A "Alex Safatli" -V ${VER} -R ${VER} -o ${DOCS} ${PKG}
	@cd ${DOCS} && make html
	cp -r ${DOCS}/_build/html/* $(GH_DOC_REPO)
	@cd $(GH_DOC_REPO) && git add * && git commit -m "Version ${VER} documentation." && git push origin gh-pages

tests: ${FILES} ${TESTS}
	cd tests && ${PY} allTests.py

clean:
	-find ${DOCS} -type f -not -name api.pdf -delete
	-find . -type f -name *.pyc -delete
	-rm -r build dist ${PKG}.egg-info

check_environment:
ifndef GH_DOC_REPO
	$(error GH_DOC_REPO is undefined)
endif
