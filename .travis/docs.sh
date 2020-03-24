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
        # Download links are generated under the assumption that the prompts in
        # the release script are followed and that the documentation artifacts
        # will be uploaded to the Github release
        echo "[pdf]: https://github.com/qucontrol/krotov/releases/download/$TRAVIS_TAG/krotov-$TRAVIS_TAG.pdf" > docs/_build/html/_downloads
        echo "[html]: https://github.com/qucontrol/krotov/releases/download/$TRAVIS_TAG/krotov-$TRAVIS_TAG.zip" >> docs/_build/html/_downloads
    else
        case "$TRAVIS_BRANCH" in
            master|doctr) DEPLOY_DIR="$TRAVIS_BRANCH";;
            *)      echo "Not deploying branch $TRAVIS_BRANCH (not in whitelist)";;
        esac
    fi
    if [ ! -z "$DEPLOY_DIR" ]; then
        python -m doctr deploy --key-path docs/doctr_deploy_key.enc \
            --command=doctr-versions-menu \
            --built-docs docs/_build/html --no-require-master --build-tags "$DEPLOY_DIR"
    fi
else
    echo "Not deploying to gh-pages (PR or not on Travis)"
fi

echo "# DOCTR - DONE"
