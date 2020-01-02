echo "# DOCTR - deploy documentation"
echo "## Generate main html documentation"
tox -c tox-pyenv.ini -e docs

if [ ! -z "$TRAVIS_TAG" ]; then

    echo "Building as tag '$TRAVIS_TAG'"

    echo "## Generate documentation downloads"
    mkdir docs/_build/download

    echo "### [htmlzip]"
    tox -c tox-pyenv.ini -e docs -- -b html _build_htmlzip
    cd docs || exit
    mv  _build_htmlzip krotov.html
    zip -r krotov.html.zip ./krotov.html
    cd ../ || exit
    mv docs/krotov.html.zip docs/_build/download

else

    echo "Building as branch '$TRAVIS_BRANCH'"

    tox -c tox-pyenv.ini -e docs -- -b latex _build/tex
    cp docs/krotovscheme.pdf docs/oct_decision_tree.pdf docs/_build/tex/
    tox -c tox-pyenv.ini -e run-cmd -- python docs/build_pdf.py
    echo "Uploading log" # DEBUG
    response=$(curl --upload-file "docs/_build/tex/krotov_final.log" "https://paste.c-net.org/") # DEBUG
    echo "$response" # DEBUG
    echo "Uploading pdf" # DEBUG
    response=$(curl --upload-file "docs/_build/tex/krotov.pdf" "https://paste.c-net.org/") # DEBUG
    echo "$response" # DEBUG
    echo "Uploaded" # DEBUG

fi


# Deploy
if [ ! -z "$TRAVIS" ]; then
    echo "## pip install doctr"
    python -m pip install doctr
    echo "## doctr deploy"
    if [ ! -z "$TRAVIS_TAG" ]; then
        DEPLOY_DIR="$TRAVIS_TAG"
    else
        case "$TRAVIS_BRANCH" in
            master) DEPLOY_DIR="$TRAVIS_BRANCH";;
            *)      echo "Not deploying branch $TRAVIS_BRANCH (not in whitelist)";;
        esac
    fi
    if [ ! -z "$DEPLOY_DIR" ]; then
        python -m doctr deploy --key-path docs/doctr_deploy_key.enc \
            --command="git show $TRAVIS_COMMIT:.travis/docs_post_process.py > post_process.py && git show $TRAVIS_COMMIT:.travis/versions.py > versions.py && python post_process.py" \
            --built-docs docs/_build --no-require-master --build-tags "$DEPLOY_DIR"
    fi
fi

echo "# DOCTR - DONE"
