import React from "react";
import type { NotificationSeverity } from "../features/notifications/types";

// Extend the Notification interface with a "level" property
export interface Notification {
  id: number;
  content: string;
  timestamp: Date;
  isNew: boolean;
  level: NotificationSeverity;
}

export interface NotificationContextType {
  notifications: Notification[];
  pushNotification: (content: string, level?: NotificationSeverity) => void;
  markAllAsRead: () => void;
  clearNotifications: () => void;
}

const NotificationContext = React.createContext<NotificationContextType>({
  notifications: [],
  pushNotification: () => {},
  markAllAsRead: () => {},
  clearNotifications: () => {},
})

export function useNotifications() {
  return React.useContext(NotificationContext);
}

export const NotificationProvider: React.FC<{ children: React.ReactNode }> = ({ children }) => {
  const [notifications, setNotifications] = React.useState<Notification[]>([]);

  const pushNotification = (content: string, level: NotificationSeverity = "info") => {
    console.log("Pushing notification:", content, level);
    setNotifications((prev) => {
      const newNotification: Notification = {
        id: Date.now(), // or use some uuid mechanism
        content,
        timestamp: new Date(),
        isNew: true,
        level, // store the level so you can adjust styling
      };

      let updated = [...prev, newNotification];

      // Keep a maximum of 50 notifications
      if (updated.length > 50) {
        updated = updated.slice(updated.length - 50);
      }
      return updated;
    })
  }

  const markAllAsRead = () => {
    setNotifications((prev) =>
      prev.map((notification) => ({ ...notification, isNew: false }))
    )
  }

  const clearNotifications = () => {
    setNotifications([]);
  }

  return (
    <NotificationContext.Provider value={{ notifications, pushNotification, markAllAsRead, clearNotifications }}>
      {children}
    </NotificationContext.Provider>
  )
}