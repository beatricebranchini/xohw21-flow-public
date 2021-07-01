# FLOW - Wavefront Alignment algorithm FPGA implementation

Team number: xohw21-211 <br />
Project name: FLOW - WFA on FPGA <br />

### Structure of the repository
src: contains source files <br />
dataset: contains an example of input file supported by our implementation <br />
original_sw: contains files of the original software version of the WFA algorithm <br />
bitstream: contains .xclbin files <br />

### Building the original software

Since we integrated the original software inside the host, it is fundamental to build the application:
```
$> cd original_sw
$> make
$> cd ..
```

### Running our implementation on FPGA

Firstly, run the following command to build the host: 
```
$> make host

```

Then, run the following commands to run everything on FPGA: <br />

1) if you want to align randomly generated inputs, run:
```
$> ./host bitstream/wfa.xclbin

```
2) if you want to align inputs given from a file, run:
```
$> ./host bitstream/wfa.xclbin datasets/sample.dataset.seq 

```
