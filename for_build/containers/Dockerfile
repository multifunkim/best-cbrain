FROM multifunkim/mcr:latest

LABEL maintainer="Obai Bin Ka'b Ali @aliobaibk"
LABEL version="1.1"


# Environment
ENV BEST_ROOT=/usr/local/Multi_FunkIm/BEst \
    BEST_DATA_DIR="$BEST_ROOT"/data \
    BEST_APP_DIR="$BEST_ROOT"/application


# Copy BEst app and test data
# (files could be hosted in a server in the future for auto build)
COPY best-data.zip /tmp/
COPY best-app.zip /tmp/


# Install BEst
RUN mkdir -p "$BEST_ROOT" && \
    mkdir "$BEST_DATA_DIR" "$BEST_APP_DIR" && \
    unzip -q /tmp/best-data.zip -d "$BEST_DATA_DIR" && \
    unzip -q /tmp/best-app.zip -d "$BEST_APP_DIR" && \
    rm -rf /tmp/* && \
    chmod +x -R "$BEST_APP_DIR"


# Update the path
ENV PATH="$BEST_APP_DIR:$PATH"