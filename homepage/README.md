# Building a Custom Homepage

If you're a fancy person who likes pretty things, you might feel the urge to spruce up the [default documentation](https://dms-vep.org/dms-vep-pipeline-3/). Below, you'll find instructions for building and deploying your own custom [VitePress](https://vitepress.dev/) website to give your documentation some [pizazz](https://dms-vep.org/Nipah_Malaysia_RBP_DMS/)! The default documentation will still be available at a different URL ([see details below](#overview)).

## Overview

You'll have to make some slight modifications to your project repo to host a [VitePress](https://vitepress.dev/) website. There are five steps:

1. [Installing the necessary packages](#installing-the-necessary-packages)
2. [Creating a `homepage/` directory](#creating-a-homepage-directory)
3. [Configuring the pipeline](#configuring-the-pipeline)
4. [Developing your site](#developing-your-site)
5. [Switching the GitHub Pages source](#switching-the-github-pages-source)

In the following sections, I'll walk you through each of these steps.

## Installing the necessary packages

You'll need to install the Javascript packages used to build and develop the [VitePress](https://vitepress.dev/) site. We'll use the 'node package manager' (`npm`) to install these packages. You can install `npm` in your existing `conda` environment by running:

```bash
conda install -c conda-forge nodejs
```

Now you'll need a `package.json` file to tell `npm` which packages to install. Copy the `package.json` file from `dms-vep-pipeline-3` into the root (top level) of your project directory. Once it's there, run the following command to install the necessary Javascript packages:

```bash
npm install
```

You should see that a `package-lock.json` file has been added to your repo.

## Creating a `homepage/` directory

Next, you'll need to add the source code for the [VitePress](https://vitepress.dev/) site. This source code will live in a new `homepage/` directory at the root of your project. Copy the example `homepage/` directory from the `dms-vep-pipeline-3` repo into your project. It should have the following structure:

```bash
homepage
├── .vitepress
│   ├── config.mjs
│   └── theme
│       ├── Altair.vue
│       ├── Figure.vue
│       ├── index.js
│       ├── parseVegaSpec.js
│       └── style.css
├── README.md
├── antibody_escape.md
├── cell_entry.md
├── index.md
├── pipeline_information.md
└── public
    ├── appendix.html
    ├── htmls
    │   ├── ...
    ├── notebooks
    │   ├── ...
    └── your-receptor.png
```

There might be some additional files – don't worry about these. The key is that you have a directory called `homepage/` at the top of your project repo with this file organization.

Now, optionally, you can test if [VitePress](https://vitepress.dev/) is working by building a site from the example code in the newly added `homepage/` directory. Make sure that you're in a `conda` environment with `npm` already installed. If you're working on the server, run the following command to boot up a local version of the example website:

```bash
npm run remote:docs:dev
```

Drop the `remote:` part of the command if you're not on the server.

```bash
npm run docs:dev
```

Now you'll see a "development server" booted up with a link that looks something like this:

```bash
http://localhost:5173/my_virus_dms/
```

This link will look different depending on whether you're running [VitePress](https://vitepress.dev/) locally or on a remote server. Either way, if you click on the link it should bring you to a local 'development' version of the website. There will be more detail on this in the [development section](#developing-your-site) of this guide. But, for now, you know that everything is working as it should!

## Configuring the pipeline

Now that you know [VitePress](https://vitepress.dev/) is running, you'll need to configure the pipeline to replace the example files in the `homepage/` directory with the output of your analysis.

### Clean up the `homepage/` directory

Remove the example files from the `homepage/` directory:

```bash
rm -rf homepage/public/htmls
rm -rf homepage/public/notebooks
rm -f homepage/public/appendix.html
```

These files will be replaced by the pipeline with the contents of your `docs/` directory.

### Modify the `.gitignore`

You need to modify the `.gitignore` to ignore the files created by [VitePress](https://vitepress.dev/). Add the following lines if they're not already present.

```bash
node_modules/
!homepage/.vitepress/
homepage/.vitepress/cache/
homepage/.vitepress/dist/
```

### Update your `config.yml`

The pipeline will automatically populate your `homepage/public` directory with the contents of your `docs/` directory. This has two benefits; it adds the default documentation as part of your new [VitePress](https://vitepress.dev/) site, and it lets you include your notebooks and `Altair` plots on the site. To tell the pipeline to do this, you'll need to update the `config.yaml' file with the following lines:

```yaml
homepage: ./homepage/public
build_vitepress_homepage: true
```

### Add a deployment workflow

You'll need to add a GitHub Actions workflow to automatically build your site when you push changes to the `main` branch of your GitHub repo. Copy the `deploy.yaml` file from `.github/workflows/deploy.yaml` in the `dms-vep-pipeline-3` repo.

### Run the pipeline

Finally, you can run the pipeline to populate the `/homepage/public` with the contents of your default documentation.

## Developing your site

## Switching the GitHub Pages source
