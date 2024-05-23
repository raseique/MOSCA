FROM continuumio/miniconda3

# the creation of a named_volume is neccessary for MP spectra conversion
VOLUME ["named_volume"]

RUN git clone https://github.com/iquasere/MOSCA.git -b development \
&& bash MOSCA/cicd/install.bash --mosca_env "$(conda info --base)" \
&& conda clean --all \
# all commands ahead are only for MP spectra conversion with docker
&& apt-get update \
&& apt-get install -y curl \
&& curl -fsSL https://get.docker.com -o get-docker.sh \
&& sh get-docker.sh

CMD [ "mosca" ]