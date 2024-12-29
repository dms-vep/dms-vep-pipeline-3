# Building a Custom Homepage

To make the results more accessible, you can improve the [default documentation](https://dms-vep.org/dms-vep-pipeline-3/) of the pipeline.
Below, you'll find instructions for building and deploying your own custom [VitePress](https://vitepress.dev/) website where you can customize the display and show key plots and data more clearly.
[Here is an example](https://dms-vep.org/Flu_H5_American-Wigeon_South-Carolina_2021-H5N1_DMS) of such a page.
The default documentation will still be available at a different URL ([see details below](#overview)).

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

Alternatively, if you're working on `Rhino` at FHCC, you can load node as a module:

```bash
ml nodejs/20.9.0-GCCcore-13.2.0
```

You'll need to make a `package.json` file to tell `npm` which packages to install. Copy the `package.json` file from `dms-vep-pipeline-3` into the root (top level) of your project directory and run the following command to install the necessary Javascript packages:

```bash
npm install
```

You should see that a `package-lock.json` file has been added to your repo.

**Note**, it's important that the major version of node that generates `package-lock.json` (with `npm install`) is same version used in [`.github/workflows/deploy.yml`](../.github/workflows/deploy.yaml) (if you don't have this file, continue to read below; you need to create it to deploy your homepage), i.e.

```yaml
strategy:
  matrix:
    node-version: [22.x]
```

Check your version of node with `node -v`. If the local version of node is different than [`.github/workflows/deploy.yml`](../.github/workflows/deploy.yaml) either upgrade your local version or upgrade the version in [`.github/workflows/deploy.yml`](../.github/workflows/deploy.yaml). It may also be necessary to run `npm audit fix` to fix your `package-lock.json`.

## Creating a `homepage/` directory

Next, you'll need to add the code for your [VitePress](https://vitepress.dev/) site. This code will live in a new `homepage/` directory at the root of your project. Copy the example `homepage/` directory from the `dms-vep-pipeline-3` repo into your project. It should have the following structure:

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
    └── your-vep.png
```

There might be some additional files – don't worry about these. The key is that you have a directory called `homepage/` at the top of your project repo with this file organization.

You can test whether [VitePress](https://vitepress.dev/) is working by building a site from the example code in the newly added `homepage/` directory. Make sure that you're in a `conda` environment with `npm` already installed. If you're working on the server, run the following command to boot up a local version of the example website:

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

The link will look different depending on whether you're running [VitePress](https://vitepress.dev/) locally or on a remote server. Either way, if you click on the link it should open a 'development' version of the website in your browser. Check out the [development section](#developing-your-site) for more details. But, for now, you know that everything is working as it should!

## Configuring the pipeline

Now that you know [VitePress](https://vitepress.dev/) is running, you'll need to configure the pipeline to replace the example files in the `homepage/` directory with the output of your analysis.

### Clean up the `homepage/` directory

Remove the example files from the `homepage/` directory:

```bash
rm -rf homepage/public/htmls
rm -rf homepage/public/notebooks
rm -f homepage/public/appendix.html
```

These files will be replaced with the contents of your `docs/` directory by the `Snakemake` pipeline.

### Modify the `.gitignore`

You need to modify the `.gitignore` to ignore certain files created by [VitePress](https://vitepress.dev/). Add the following lines to your `.gitignore` if they're not already present.

```bash
node_modules/
!homepage/.vitepress/
homepage/.vitepress/cache/
homepage/.vitepress/dist/
```

### Update your `config.yml`

The pipeline will automatically populate your `homepage/public` directory with the contents of your `docs/` directory. This has two benefits; it adds the default documentation as part of your new [VitePress](https://vitepress.dev/) site, and it lets you include your notebooks and `Altair` plots on the site. To tell the pipeline to do this, you'll need to update the pipeline's `config.yaml` file with the following lines:

```yaml
homepage: homepage/public
build_vitepress_homepage: true
```

### Add a deployment workflow

You'll need to add a GitHub Actions workflow to automatically build your site when you push changes to the `main` branch of your GitHub repo. Copy the file [`.github/workflows/deploy.yaml`](../.github/workflows/deploy.yaml) in the `dms-vep-pipeline-3` repo into your `.github/workflows` directory in your repository (you may have to create this directory).

### Run the pipeline

Finally, you can run the pipeline to populate the `homepage/public` with the contents of your default documentation.

## Developing your site

At this point, you can start adding content to your fancy new [VitePress](https://vitepress.dev/) site! Don't worry about making your unfinished content public, the new site won't be visible on GitHub Pages until the final step.

Making a [VitePress](https://vitepress.dev/) site is as simple as writing Markdown documents. Any file ending in `.md` in the `homepage/` directory is interpreted by [VitePress](https://vitepress.dev/) as a page on your website. All of the features of Markdown are available alongside features enabled by `html` and custom 'components' discussed in detail [below](#custom-components).

### Configuring your site

You'll need to replace the default content in the [VitePress](https://vitepress.dev/) config file (`homepage/.vitepress/config.mjs`).

```js
// https://vitepress.dev/reference/site-config
export default defineConfig({
  lang: "en-US",
  title: "Homepage for my favorite viral entry protein",
  description:
    "Data, figures, and analysis for DMS of my favorite viral entry protein.",
  base: "/my_virus_dms/",
  themeConfig: {
    nav: [
      { text: "Home", link: "/" },
      { text: "Appendix", link: "/appendix", target: "_self" },
    ],
    socialLinks: [{ icon: "github", link: "https://github.com/dms-vep" }],
    footer: {
      message: "Copyright © 2024-present Me and Jesse Bloom",
    },
  },
});
```

Replace the `title` and `description` with content that suits your project. Add your name to the copyright message in the `footer`. Change the `link` under `socialLinks` to your project's GitHub repo. Finally, replace the `base` path with the name of your repository. Setting the `base` path correctly is essential for hosting your site.

### Configuring the landing page

The landing page for your [VitePress](https://vitepress.dev/) site is determined by the `index.md` file. Edit this file to change the title add an image of your viral entry protein and link to pages describing your results.

```yml
layout: home

hero:
  name: "My Favorite Virus' Receptor DMS"
  tagline: "A collection of data, figures, and information for DMS of my favorite viral entry protein"
  image: your-vep.png
features:
  - title: Antibody Escape
    details: Example link to antibody escape data
    link: /antibody_escape
  - title: Cell Entry
    details: Example link to functional selection data
    link: /cell_entry
  - title: Pipeline Information
    details: Example page detailing computational pipeline
    link: /pipeline_information
```

The `index.md` file, like every Markdown page of your site, is comprised of two parts: optional `yaml` 'frontmatter' and Markdown content. In the case of the `index.md`, the content of the page is dictated solely by the `yaml` 'frontmatter'. In the rest of the pages, you'll write Markdown content with some minor configuration through the `yaml` 'frontmatter'. 

Check out the [VitePress Writing guide](https://vitepress.dev/guide/markdown) for details on writing pages and configuring them with 'frontmatter'.

### Writing pages

You add pages to your website by writing Markdown documents in the `homepage/` directory. VitePress will convert these Markdown documents into `html` and render them as pages on your website. You can link to these pages in your `index.md` (or anywhere else on the site) using their path with the `.md` extension omitted. For example, if you write a page with details about antibody escape called `antibody_escape.md`, the `/antibody_escape` path will bring you to that page.

### Key syntax

As mentioned above, each page is a simple [Markdown document](https://vitepress.dev/guide/markdown) configured with[ `yaml` frontmatter](https://vitepress.dev/guide/frontmatter). However, there are some key points to keep in mind.

For one, if you reference a link to an `html` file created by the pipeline like an `Altair` plot, you'll need to tell VitePress that you're linking to a *local* file and not a *remote* URL. You do this by adding the `target="_self"` directive after the relevant Markdown link.

```md
[LibA-220210-REGN10933-1](notebooks/fit_escape_antibody_escape_LibA-220210-REGN10933-1.html){target="_self"}
```

You can also render interactive `Altair` plots directly on a page with a custom component that I added to `homepage/.vitepress`.

```html
<Figure caption="Effects of mutations on antibody neutralization by REGN10933">
    <Altair :showShadow="true" :spec-url="'htmls/REGN10933_mut_effect.html'"></Altair>
</Figure>
```

The code above will render the `Altair` plot `htmls/REGN10933_mut_effect.html` in an expandable widget on the page with the figure caption "Effects of mutations on antibody neutralization by REGN10933" beneath it.

## Switching the GitHub Pages source

Now that you're happy with the content of your new [VitePress](https://vitepress.dev/) site, you can make it public and show the world. Currently, GitHub Pages is still hosting the original documentation located in the `docs/` directory. However, your [VitePress](https://vitepress.dev/) site is being built in a `gh-pages` branch by GitHub Actions each time you push to the `main` branch. By switching the GitHub Pages *source* from the `docs/` directory to the `gh-pages` branch, GitHub Pages will begin serving the [VitePress](https://vitepress.dev/)site instead of the original documentation. But don't worry! Your original documentation will still be available on the new [VitePress](https://vitepress.dev/) site in a tab on the top of the page called `appendix`.

## Adding Google Analytics
If you want to add Google Analytics to your site, you can do this to see how many people have visited it.
To do this, first you need to use your Google account to create a Google Analytics tag.
Once you have done that, you can add it to the `.vitepress/config.mts` file as indicated below (in that example, `G-P7HL8Q4F41` is the Google Analytics tag; you will get a different tag when you create yours:

      head: [
        [
          "script",
          { async: "", src: "https://www.googletagmanager.com/gtag/js?id=G-P7HL8Q4F41" }
        ],
        [
          "script",
          {},
          `window.dataLayer = window.dataLayer || [];
          function gtag(){dataLayer.push(arguments);}
          gtag("js", new Date());
          gtag("config", "G-P7HL8Q4F41");`
        ]
      ],
  themeConfig: {
