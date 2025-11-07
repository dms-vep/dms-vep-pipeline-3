import { defineConfig } from "vitepress";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  lang: "en-US",
  title: "Homepage for `dms-vep-pipeline-3` test example",
  description:
    "Data, figures, and analysis for `dms-vep-pipeline-3` test example.",
  base: "/dms-vep-pipeline-3/",
  appearance: false,
  themeConfig: {
    // https://vitepress.dev/reference/default-theme-config
    nav: [
      { text: "Home", link: "/" },
      { text: "Appendix", link: "/appendix", target: "_self" },
    ],
    socialLinks: [{ icon: "github", link: "https://github.com/dms-vep" }],
    footer: {
      message: "Copyright © 2024-present Will Hannon and Jesse Bloom",
    },
  },
});
