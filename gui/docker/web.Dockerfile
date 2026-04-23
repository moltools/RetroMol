# Build the React app and serve with nginx
FROM node:20-alpine AS build
WORKDIR /usr/src/app

# Client code lives in src/client
COPY gui/src/client/package*.json ./
RUN npm ci
COPY gui/src/client ./
RUN npm run build

# nginx stage
FROM nginx:alpine
# nginx config (proxies /api to backend)
COPY gui/docker/nginx.conf /etc/nginx/conf.d/default.conf
# Static build
COPY --from=build /usr/src/app/build /usr/share/nginx/html
EXPOSE 80