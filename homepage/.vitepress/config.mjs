import { defineConfig } from "vitepress";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  lang: "en-US",
  title: "Homepage for my favorite viral entry protein",
  description:
    "Data, figures, and analysis for DMS of my favorite viral entry protein.",
  base: "/my_virus_dms/",
  themeConfig: {
    // https://vitepress.dev/reference/default-theme-config
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
