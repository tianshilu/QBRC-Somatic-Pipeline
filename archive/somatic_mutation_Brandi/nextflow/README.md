# Astrocyte Example Workflow Package

This is an example workflow package for the BioHPC Astrocyte workflow engine. Astrocyte is a system allowing workflows to be run easily from the web, in a push-button manner, taking advantage of the BioHPC compute cluster. Astrocyte allows users to access this workflow package using a simple web interface, created automatically from the definitions in this package.

This workflow package provides:

  1)  A sample ChIP-Seq data analysis workflow, which uses BWA and MACS to call peaks from one or more ChIP-Seq FASTQ input files. The workflow is implemented in the *Nextflow* workflow language.
  
  2) A sample *R Shiny* visualization app, which provides a web-based tool for visualizing results.
  
  3) Meta-data describing the workflow, it's inputs, output etc.
  
  4) User-focused documentation, in markdown format, that will be displayed to users in the Astrocyte web interface.
  
  5) Developer-focused documentation, in this file.
  
## Workflow Package Layout

Workflow packages for Astrocyte are Git repositories, and have a common layout which must be followed so that Astrocyte understands how to present them to users.


### Meta-Data

  * `astrocyte_pkg.yml` - A file which contains the metadata describing the workflow in human & machine readable text format called *YAML*. This includes information about the workflow package such as it's name, synopsis, input parameters, outputs etc.


### The Workflow
  
  * `workflow/main.nf` - A *Nextflow* workflow file, which will be run by Astrocyte using parameters provided by the user.
  * `workflow/scripts` - A directory for any scripts (e.g. bash, python, ruby scripts) that the `main.nf` workflow will call. This might be empty if the workflow is implemented entirely in nextflow. You should *not* include large pieces of software here. Workflows should be designed to use *modules* available on the BioHPC cluster. The modules a workflow needs will be defined in the `astrocyte_pkg.yml` metadata file.
  * `workflow/lib` - A directory for any netflow/groovy libraries that might be included by workflows using advanced features. Usualy empty.

### The Visualization App *(Optional)*

  * `vizapp/` - A directory that will contain an *R Shiny* visualization app, if required. The vizualization app will be made available to the user via the Astrocyte web interface. At minimum the directory requires the standard Shiny `ui.R` and `server.R` files. The exact Shiny app structure is not prescribed. Any R packages required by the Shiny app will be listed in the `astrocyte_pkg.yml` metadata.

  The visualization app will have access to any final output that was published to the `$baseDir` location in the
  nextflow workflow. This path will be accessible as `Sys.getenv('baseDir')`.


### User Documentation 

  * `docs/index.md` - The first page of user documentation, in *markdown* format. Astrocyte will display this documentation to users of the workflow package.

  * `docs/...` - Any other documentation files. *Markdown* `.md` files will be rendered for display on the web. Any images used in the documentation should also be placed here.

  
### Developer Documentation

  * `README.md` - Documentation for developers of the workflow giving a brief overview and any important notes that are not for workflow users.
  * `LICENSE.md` *(Optional)* - The license applied to the workflow package.
  * `CHANGES.md` - A brief summary of changes made through time to the workflow.

### Testing

  * `test_data/` - Every workflow package should include a minimal set of test data that allows the workflow to be run, testing its features. The `test_data/` directory is a location for test data files. Test data should be kept as small as possible. If large datasets (over 20MB total) are unavoidable provide a `fetch_test_data.sh` bash script which obtains the data from an external source.
  * `test_data/fetch_test_data.sh` *Optional* - A bash script that fetches large test data from an external source, placing it into the `test_data/` directory.

  
## Testing/Running the Workflow with the Astrocyte CLI

Workflows will usually be run from the Astrocyte web interface, by importing the workflow package repository and making it available to users. During development you can use the Astrocyte CLI scripts to check, test, and run your workflow against non-test data.

To check the structure and syntax of the workflow package in the directory `astrocyte_example`:

```bash
$ astrocyte_cli check astrocyte_example
```

To launch the workflows defined tests, against included test data:

```bash
$ astrocyte_cli test astrocyte_example
```

To run the workflow using specific data and parameters. A working directory will be created.

```bash
$ astrocyte_cli run astrocyte_example --parameter1 "value1" --parameter2 "value2"...
```

To run the Shiny vizualization app against test_data

```bash
$ astrocyte_cli shinytest astrocyte_example
```

To run the Shiny vizualization app against output from `astrocyte_cli run`, which will be in the work directory created by `run`:

```bash
$ astrocyte_cli shiny astrocyte_example
```

To generate the user-facing documentation for the workflow and display it in a web browser"
```bash
$ astrocyte_cli docs astrocyte_example
```

## The Example ChIP-Seq Workflow ##



## The Example Shiny Viz App ##

## Provided Test Data ##