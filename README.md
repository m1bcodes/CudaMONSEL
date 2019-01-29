# CudaMONSEL
Monte Carlo Simulation of SEM Signals using CUDA.

## Citing CudaMONSEL:

If you use CudaMONSEL in your research, please cite with:
```
@misc{villarrubia2015jmonsel,
  title={Scanning electron microscope measurement of width and shape of 10 nm patterned lines using a JMONSEL-modeled library},
  author={J.S. Villarrubia et. al.},
  howpublished={\url{https://ws680.nist.gov/publication/get_pdf.cfm?pub_id=916512}},
  year={2015}
}
```

## Common Problems
### "unspecified launch failure"
https://stackoverflow.com/questions/36903282/why-does-my-cuda-kernel-crash-unspecified-launch-failure-with-a-different-data

### "CUDALINK : nvlink error : Undefined reference to '*' in 'Debug/*.obj'"
Every function in a class must be defined for the class to be used on the device.

https://stackoverflow.com/questions/37507274/getting-error-nvlink-error-undefined-reference-to-zn8strategy8backtestepdd
