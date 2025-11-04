echo "construyendo programas"
python constructor.py
echo "programas construidos"

echo "compilando librerias"
gfortran -c metaheuristicas.f90
gfortran -c tsukamoto.f90
gfortran -c OptimizationDriver.f90
gfortran metaheuristicas.o tsukamoto.o OptimizationDriver.o -o driver.exe
echo "librerias compiladas"

echo "ejecutando driver"
driver.exe
echo "simulación terminada"

echo "Comienza verificaciones"
python verificaciones.py

gfortran -c verificaciones.f90
gfortran tsukamoto.o verificaciones.o -o verificaciones.exe
verificaciones.exe

python unir_resultados.py

python datos_vv.py

gfortran -c verificaciones.f90
gfortran tsukamoto.o verificaciones.o -o verificaciones.exe
verificaciones.exe

echo "Realizando limpieza"
del datos_vv.py
del driver.exe
del metaheuristicas.f90
del metaheuristicas.mod
del metaheuristicas.o
del OptimizationDriver.f90
del OptimizationDriver.o
del rs.csv
del tsukamoto.f90
del tsukamoto.mod
del tsukamoto.o
del unir_resultados.py
del verificaciones.exe
del verificaciones.f90
del verificaciones.o
del verificaciones.py
