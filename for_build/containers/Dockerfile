FROM multifunkim/matlab-runtime:r2019b

LABEL maintainer="@aliobaibk" \
      version="1.1" \
      description="BEst standalone version"



# Environment
ENV LANG='C.UTF-8' \
    MFIL_HOME="/usr/local/Multi_FunkIm" \
    BEST_VERSION="cbrain-mcv97"
ENV BEST_DIR="$MFIL_HOME/best_$BEST_VERSION"
ENV BEST_DATA_DIR="$BEST_DIR/data"



# Install BEst (files could be hosted in a server in the future for auto build)
COPY best-bin.zip /tmp/
COPY best-data.zip /tmp/
RUN mkdir -p "$MFIL_HOME/bin" "$BEST_DIR" && \
    unzip -q /tmp/best-bin.zip -d "$BEST_DIR" && \
    unzip -q /tmp/best-data.zip -d "$BEST_DIR/data" && \
    rm -f /tmp/best-bin.zip /tmp/best-data.zip && \
    chmod +x -R "$BEST_DIR" && \
    ln -s "$BEST_DIR/best_samapp" "$MFIL_HOME/bin/best"



# Update environment
ENV PATH="$MFIL_HOME/bin:$PATH"