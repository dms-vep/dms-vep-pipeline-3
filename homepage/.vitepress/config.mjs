import { defineConfig } from "vitepress";

// https://vitepress.dev/reference/site-config
export default defineConfig({
  lang: "en-US",
  title: "Nipah RBP DMS",
  title: "My Favorite Virus' Receptor DMS",
  description:
    "Data, figures, and analysis for DMS of my favorite viral receptor.",
  base: "/my_virus_dms/",
  themeConfig: {
    // https://vitepress.dev/reference/default-theme-config
    nav: [
      { text: "Home", link: "/" },
      { text: "Cell Entry", link: "/cell_entry" },
      { text: "Antibody Escape", link: "/antibody_escape" },
      { text: "Pipeline Information", link: "/pipeline_information" },
    ],

    socialLinks: [{ icon: "github", link: "https://github.com/dms-vep" }],
  },
});
