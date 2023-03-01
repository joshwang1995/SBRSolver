<a name="readme-top"></a>
# SBR Ray Tracer

<!-- ABOUT THE PROJECT -->
## About The Project

Implementation of the Shoot Bounce Relfect (SBR) method for electromagnetics simulation.

Features
* Transmitter Modelling
  - Ray generation via icosahedron sphere tessellation
  - Choose and apply theoretical antenna patterns (Hertzian Dipole, half-wave dipole, isotropic)
* Receiver Modelling
  - Reception sphere algorithm
  - Received field computation
* Scene Modelling
  - STL to triangle primitives
  - Bounding Volume Hierachy (BVH) implementation
  - Coplanar surface identification
* Ray Path Tracing
  - Ray math library
  - Ray launching algorithm
  - Ray intersection tests
  - Duplicate path removal
  - Image method path correction
* Field Calculation
  - TE/TM field decomposition
  - TE/TM field reflection
  - TE/TM field transmission

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- GETTING STARTED -->
## Getting Started
Follw the instructions below to download the required software and compile the source code. 

### Prerequisites
Required Software
* Visual Studio (latest version)
  - Download link: https://visualstudio.microsoft.com/
* Your GitHub account needs to be added as a Collaborator on the sarrisgroup/raytracer repository (ask Prof. Sarris to add you).

### Installation
1. Open Visual Studio &rarr; <kbd>File | Clone Repository </kbd>
2. Under <kbd>Browse a repository</kbd>, click on <kbd>GitHub</kbd>
3. If you havn't signed into your GitHub account, you can do so by clicking the <kbd>Sign In</kbd> box.
4. Under <kbd>Collaborator repositories</kbd>, choose <kbd>sarrisgroup/raytracer</kbd>
5. Pick the path to clone the repo to in the <kbd>Local Path:</kbd> box
6. Click <kbd>Clone</kbd> and project is now cloned locally
7. **Before making any changes to the code, read the Usage section**

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- USAGE EXAMPLES -->
## Usage

**_Git is used for version control in this entire proejct, any use of Git without fully understanding the consequence of each action can be detrimental. Read the following resource to get started with using Git. https://git-scm.com/book/en/v2_**

1. Once you have a local repo, you need to make sure you are not changing code that is in the master branch. This branch should only contain verified working code.
2. To make a branch from master, click <kbd>Git | New Local Branch From</kbd>
3. Choose a branch name and click <kbd>Create</kbd>
4. Now you should have your own branch where the codebase is the exact replica of the master branch. Any changes in the local repo will be recorded. 
5. It is a good idea to frequently commit your code so you can go back if something goes wrong. Think of it as quick save in video games.
  - To commit your code, click <kbd>Git | Commit or Stash...</kbd>
  - In the side panel window, write a short and concise commit message (insert link here for good commit message practices), then click <kbd>Commit All</kbd>
6. You can push the commit to the remote branch by clicking the <kbd> Git | Push </kbd>. This uploades your code changes from your local repo to the GitHub repo
7. Once you finish your code fixes/implementation on your remote branch, you can create a pull request to merge your code changes to the master branch. **Make sure you test out your code as much as you can before merging and get someone else to review it for you**

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- TROUBLESHOOTING -->
## Troubleshooting

* When you clone the entire project from Visual Studio, it should recursively clone all submodules. If the submodule folder is empty, run the following command in the directory where the git repo exists
   ```sh
   git submodule init & git submodule update
   ```
<p align="right">(<a href="#readme-top">back to top</a>)</p>


