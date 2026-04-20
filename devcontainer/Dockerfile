FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN set -eux; \
    sed -i -e 's|http://archive.ubuntu.com/ubuntu|https://archive.ubuntu.com/ubuntu|g' \
           -e 's|http://security.ubuntu.com/ubuntu|https://security.ubuntu.com/ubuntu|g' /etc/apt/sources.list; \
    printf 'Acquire::Retries "5";\nAcquire::http::Timeout "30";\nAcquire::https::Timeout "30";\n' > /etc/apt/apt.conf.d/99retries

# --- system deps ---
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget gcc-14 g++-14 \
    sudo \
    libgl1 libglib2.0-0 libsm6 libxrender1 libxext6 libegl1 \
    libx11-xcb1 libxcb-cursor-dev libxcb-icccm4 libxcb-image0 libxcb-keysyms1 \
    libxcb-render-util0 libxcb-render0 libxcb-shape0 libxcb-xfixes0 libxcb-sync1 \
    libxcb-util1 libxcb-xinerama0 libxrandr2 libdbus-1-3 \
    xauth x11-utils x11-xserver-utils qt6-base-dev xvfb \
    libz-dev xserver-xorg-core zlib1g-dev build-essential libc6-dev libstdc++6 \
    jellyfish \
 && rm -rf /var/lib/apt/lists/*

# --- MS core fonts (optional) ---
RUN wget http://ftp.br.debian.org/debian/pool/contrib/m/msttcorefonts/ttf-mscorefonts-installer_3.8_all.deb && \
    apt update && apt install -y --no-install-recommends fontconfig ./ttf-mscorefonts-installer_3.8_all.deb && \
    rm -f ttf-mscorefonts-installer_3.8_all.deb && fc-cache -f -v && \
    rm -rf /var/lib/apt/lists/*

# --- Miniconda ---
RUN wget -q https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /miniconda.sh && \
    bash /miniconda.sh -b -p /opt/conda && rm /miniconda.sh
ENV PATH="/opt/conda/bin:${PATH}"
SHELL ["/bin/bash", "-c"]

# conda初期化（グローバル）
RUN . /opt/conda/etc/profile.d/conda.sh && conda init bash

# ★ "conda tos accept" は存在しないので削除 ★

# --- mamba + env ---
COPY environment.yml /tmp/environment.yml
RUN conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r
RUN conda install -y -n base -c conda-forge mamba && \
    mamba env create -f /tmp/environment.yml && \
    conda clean --all -y && rm -f /tmp/environment.yml

# すべてのユーザーで dev を自動有効化
RUN echo '. /opt/conda/etc/profile.d/conda.sh' >> /etc/bash.bashrc && \
    echo 'conda activate dev' >> /etc/bash.bashrc

ENV QT_DEBUG_PLUGINS=1

# --- Set up a non-root user (修正版) ---
ARG USERNAME=saori
ARG USER_UID=1000
ARG USER_GID=1000

RUN set -eux; \
    if id -u ubuntu >/dev/null 2>&1; then userdel -r ubuntu || true; fi; \
    # 既存GIDがあればそのグループ名をUSERNAMEに変更、なければ作成
    if getent group "${USER_GID}" >/dev/null; then \
        groupmod -n "${USERNAME}" "$(getent group "${USER_GID}" | cut -d: -f1)"; \
    else \
        groupadd --gid "${USER_GID}" "${USERNAME}"; \
    fi; \
    # ユーザー作成/調整
    if id -u "${USERNAME}" >/dev/null 2>&1; then \
        usermod -u "${USER_UID}" -g "${USER_GID}" -s /bin/bash -d "/home/${USERNAME}" "${USERNAME}"; \
    else \
        useradd --uid "${USER_UID}" --gid "${USER_GID}" -m -s /bin/bash "${USERNAME}"; \
    fi; \
    usermod -aG sudo "${USERNAME}"; \
    install -d -m 0755 /etc/sudoers.d; \
    echo "${USERNAME} ALL=(root) NOPASSWD:ALL" | tee "/etc/sudoers.d/${USERNAME}" >/dev/null; \
    chmod 0440 "/etc/sudoers.d/${USERNAME}"; \
    echo '. /opt/conda/etc/profile.d/conda.sh' >> "/home/${USERNAME}/.bashrc"; \
    echo 'conda activate dev' >> "/home/${USERNAME}/.bashrc"; \
    mkdir -p /home/result/tmp; \
    chown -R "${USERNAME}:${USERNAME}" /home/${USERNAME} /home/result

ENV ETERNAFOLD_PATH=/opt/conda/envs/dev/bin/eternafold-bin \
    ETERNAFOLD_PARAMETERS=/opt/conda/envs/dev/lib/eternafold-lib/parameters/EternaFoldParams.v1 \
    VIENNA_2_PATH=/opt/conda/envs/dev/bin \
    CONTRAFOLD_PATH=/opt/conda/envs/dev/bin/eternafold-bin \
    LINEARFOLD_PATH=/opt/conda/envs/dev/bin \
    TMP=/home/result/tmp

USER ${USERNAME}
CMD ["bash"]
