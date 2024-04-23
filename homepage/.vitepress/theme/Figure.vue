<script>
export default {
    props: {
        imageSrc: {
            type: String,
            default: ''
        },
        caption: {
            type: String,
            required: true
        }
    },
};
</script>

<template>
    <p>
    <figure>
        <img v-if="imageSrc" :src="imageSrc" :alt="caption" />
        <slot v-else></slot>
        <figcaption>{{ caption }}</figcaption>
    </figure>
    </p>
</template>
  
<style scoped>
figcaption {
    color: grey;
    padding-top: 0.5em;
    font-size: 0.8em;
    text-align: center;
}
</style>