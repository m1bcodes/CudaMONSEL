# CudaMONSEL
<!--
![alt text](https://raw.githubusercontent.com/zengrz/CudaMONSEL/double-to-float/CudaMONSEL/12.png)
-->
<p align="center"> 
<img src="https://raw.githubusercontent.com/zengrz/CudaMONSEL/double-to-float/CudaMONSEL/12.png"> 
</p>

CudaMONSEL is a full-fledged electron tracker based on first physical principles. Its primary application is to carry out Monte Carlo simulation of SEM Signals. It can be ran on CPU using a thread pool, as well as GPU using the CUDA framework.

CudaMONSEL is a direct port of JMONSEL, a Java version the software built by J.S. Villarrubia and Nicholas Ritchie of NIST. CudaMONSEL aims to speed up the original simulation to mass produce SEM images for ML training. It also added extra functionalities to describe the geometry of the setup. For example, 2D projection of outline of shapes on an arbitrary plane.

## Citing CudaMONSEL:

If you use CudaMONSEL in your research, please cite with:
```
@misc{villarrubia2015jmonsel,
  title={Scanning electron microscope measurement of width and shape of 10 nm patterned lines using a JMONSEL-modeled library},
  author={J.S. Villarrubia et. al.},
  howpublished={\url{https://ws680.nist.gov/publication/get_pdf.cfm?pub_id=916512}},
  year={2015}
}

@misc{zeng2019cudamonsel,
  title={CudaMONSEL},
  author={Ruizi, Zeng},
  howpublished={\url{https://github.com/zengrz/CudaMONSEL/}},
  year={2019}
}
```
## Citing Amphibian:
Amphibian is an initiative to build a library of data structures and algorithms that can be used on both host and device. Can be useful when transitioning from CPU to CUDA code.

```
@misc{zeng2019amphibian,
  title={Amphibian},
  author={Zeng, Ruizi},
  howpublished={\url{https://github.com/zengrz/Amphibian/}},
  year={2019}
}
```
