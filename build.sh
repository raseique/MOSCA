echo "Storing MOSCA's files in the Conda environment at: ${PREFIX}"
mkdir -p "${PREFIX}/share/MOSCA" "${PREFIX}/bin"
cp -r MOSCA/workflow/* "${PREFIX}/share/MOSCA"
chmod +x "${PREFIX}/share/MOSCA/mosca.py"
ln -s "${PREFIX}/share/MOSCA/mosca.py" "${PREFIX}/bin/mosca"