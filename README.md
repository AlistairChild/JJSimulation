# JJSimulator

## Introduction

JJSimulator is a mathematical tool used to model the critical current through a Josephson Junction when a magnetic field is applied. A 3D model can be generated showing how the critical current varies with the direction and strength of the applied field. Changes to the junctions charge density distributions can be modelled allowing possible insight into certain junctions charge distributions (also creates a great looking graph!)

 The mathematics is based on the mathematical theory as presented in Antonio Barone and GianFranco Paterno book "Physics And Applications Of The Josephson Effect". 

 The data generates can be output to "tsv" or "csv" and plotted accordingly. A great plotting software to use for the plotting is [plot3d](https://oomzay.github.io/plot3d/).

![Overview Animation](images/1.png)

## Usage

After creating the exicutable options can be passed in as parameters. These consist of

| Arguments      | Description |
| ----------- | ----------- |
| **-help**      | :     shows options  |
| **--lowfield**  | :     lowest field used (default 0)     |
| **--highfield**   | :     highest field used (default 0.15)      |
| **--magres**   | :     n points per angle (default 100)     |
| **--intres**   | :     n trapezoid integral! (default 100)      |
| **--length**   | :     length default (500E-9)     |
| **--width**   | :    width default (500E-9)      |
| **--height**   | :     height default (1E-9)      |
| **--angleincrement**   | :    how much to inrement angle by 0-90 (default 1)     |
| **--step_height**   | :     fraction of max critica current. (default 0.1)      |
| **--write_to**   | :     csv or tsv      |
| **--distribution**   | :     current density 1d or 2d     |
| **--view**   | :     profile or critical_current     |
| **--step_thickness**   | :     fraction of length.(default 0.1)      |

To view the current density distribution user view "profile". To view the critical current 3d plot view "critical_current".

## License

Licenced under the [MIT Licence](LICENCE.md)

Copyright (c) 2021 Alistair Child

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## Authors and acknowledgment

This project was a collaboration between Alistair Child and under a project supervised by Niladri Banerjee. 