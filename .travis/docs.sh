echo "# DOCTR - deploy documentation"
echo "## Generate main html documentation"
tox -c tox-pyenv.ini -e docs

if [ ! -z "$TRAVIS_TAG" ]; then
    echo "Building as tag '$TRAVIS_TAG'"
elif [ ! -z "$TRAVIS_BRANCH" ]; then
    echo "Building as branch '$TRAVIS_BRANCH'"
else
    echo "At least one of TRAVIS_TAG and TRAVIS_BRANCH must be set"
    sync
    exit 1
fi


# Deploy
if [ ! -z "$TRAVIS" ] && [ "$TRAVIS_EVENT_TYPE" != "pull_request" ]; then
    echo "## pip install doctr"
    python -m pip install doctr-versions-menu
    echo "## doctr deploy"
    if [ ! -z "$TRAVIS_TAG" ]; then
        DEPLOY_DIR="$TRAVIS_TAG"
    else
        case "$TRAVIS_BRANCH" in
            master|doctr) DEPLOY_DIR="$TRAVIS_BRANCH";;
            *)      echo "Not deploying branch $TRAVIS_BRANCH (not in whitelist)";;
        esac
    fi
    if [ ! -z "$DEPLOY_DIR" ]; then
        python -m doctr deploy --key-path docs/doctr_deploy_key.enc \
            --command=doctr-versions-menu \
            --built-docs docs/_build --no-require-master --build-tags "$DEPLOY_DIR"
    fi
else
    echo "Not deploying to gh-pages (PR or not on Travis)"
fi

echo "# DOCTR - DONE"
