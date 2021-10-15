<div id="top"></div>

[![MIT License][license-shield]][license-url]

<h3 align="center">Making the "Lobster" Plot</h3>

  <p align="center">
    This code contains all the ingredients to make the so-called "Lobster plot", which is the dependence of the parameter m𝛽𝛽 on the presently unknown lightest
neutrino mass, with 3-sigma bands showing allowed regions given our knowledge about the neutrino mixing parameters. This code was originally written by fromer Yale grad student Jeremy Cushman, and I updated and documented it, as well as added a couple of other functions.
    <br />
    <a href="https://github.com/toej93/LobsterPlot"><strong>Explore the docs »</strong></a>
    <br />
    <br />

  </p>
<!-- </div> -->

![Example plot](https://raw.githubusercontent.com/toej93/LobsterPlot/main/Lobster_example.png)


<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>

<!-- GETTING STARTED -->
## Getting Started

This is an example of how you may give instructions on setting up your project locally.
To get a local copy up and running follow these simple example steps.

### Prerequisites

You need the following software to be able to run this
* Python 3
* Jupyter (if you want to use the notebook, otherwise not needed)

### Installation

1. No installation needed, just clone and run!


<!-- USAGE EXAMPLES -->
## Usage

The functions needed to make plots are contained in the `nu_mass.py` file, along with their docstrings. Examples of some of the different plots that can be made are contained in [`lobsterPlot.ipynb`](https://github.com/toej93/LobsterPlot/blob/main/lobsterPlot.ipynb).


### Code to sample parameter space

I am already providing the needed sampled arrays to produce the final plots, but you can produce your own files by running `getArrays.py`. It produces two files: `normalArray.npy` and `invertedArray.npy`, which contain the sampled values as well as the 3-sigma bands.

### Sources for the included experimental limits:

* Ge limit from GERDA (2020): https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.252502
* Mo limit from CUPID-Mo (2021): https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.126.1
* Te limit from CUORE (PRELIMINARY, 2021): https://arxiv.org/abs/2104.06906
* Xe limit from KamLand-Zen (2016): https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.117.082503

<!-- CONTRIBUTING -->
## Contributing

This is an open source code, so please feel free to use it, and modify it, at your convenience. However, any contributions to this code **greatly appreciated**. 

If you have a suggestion or correction, please fork the repo and create a pull request. You can also simply open an issue with the tag "enhancement".

1. Fork the Project
2. Create your Feature Branch (`git checkout -b feature/Feature`)
3. Commit your Changes (`git commit -m 'Add some Feature'`)
4. Push to the Branch (`git push origin feature/Feature`)
5. Open a Pull Request



<!-- LICENSE -->
## License

Distributed under the MIT License. See `LICENSE.txt` for more information.




<!-- CONTACT -->
## Contact

Jorge Torres - [@toej93](https://twitter.com/toej93) - jorge.torresespinosa@yale.edu

Project Link: [https://github.com/toej93/LobsterPlot](https://github.com/toej93/LobsterPlot)


<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

* Thanks a lot to Jeremy Cushman for writing the bulk of the code.
* Thanks to the CUORE@Yale folks for reviewing the code.

<p align="right">(<a href="#top">back to top</a>)</p>



<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->

[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=for-the-badge
[license-url]: https://github.com/toej93/LobsterPlot/blob/main/LICENSE
