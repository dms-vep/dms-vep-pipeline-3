---
aside: false
---

# Antibody Escape

Here's what a page with antibody escape data might look like.

## Individual Antibody Selections

You can include links to the antibody escape data created automatically by the [`dms-vep-pipeline`](https://github.com/dms-vep/dms-vep-pipeline-3).

- [LibA-220210-REGN10933-1](notebooks/fit_escape_antibody_escape_LibA-220210-REGN10933-1.html){target="_self"}
- [LibA-220210-REGN10933-2](notebooks/fit_escape_antibody_escape_LibA-220210-REGN10933-2.html){target="_self"}
- [LibB-220302-REGN10933-1](notebooks/fit_escape_antibody_escape_LibB-220302-REGN10933-1.html){target="_self"}
- [LibA-220302-S2M11-1](notebooks/fit_escape_antibody_escape_LibA-220302-S2M11-1.html){target="_self"}
- [LibA-220302-S2M11-2](notebooks/fit_escape_antibody_escape_LibA-220302-S2M11-2.html){target="_self"}

You do this using normal markdown syntax with one exception – you include `{target="_self"}` after the link:

```bash
- [LibA-220210-REGN10933-1](notebooks/fit_escape_antibody_escape_LibA-220210-REGN10933-1.html){target="_self"}
```

This approach will work for any file that's ordinarily present in the `docs/` directory created by the pipeline. For example, `.html` files that contain `Altair` plots.

- [REGN10933 Escape](htmls/REGN10933_mut_effect.html){target="_self"}
- [S2M11 Escape](htmls/S2M11_mut_effect.html){target="_self"}

If you prefer, you can include `Altair` plots directly on the page.

<Figure caption="Effects of mutations on antibody neutralization by REGN10933">
    <Altair :showShadow="true" :spec-url="'htmls/REGN10933_mut_effect.html'"></Altair>
</Figure>

You do this using a custom html-like syntax.

```bash
<Figure caption="Effects of mutations on antibody neutralization by REGN10933">
    <Altair :showShadow="true" :spec-url="'htmls/REGN10933_mut_effect.html'"></Altair>
</Figure>
```

The code above creates the `Altair` plot using the `<Altair></Altair>` tag and adds a figure caption with the surrounding `<Figure></Figure>` tag.
